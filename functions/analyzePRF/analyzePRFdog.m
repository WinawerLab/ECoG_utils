function results = analyzePRFdog(stimulus,data,tr,options)

% function results = analyzePRFdog(stimulus,data,tr,options)
%
% <stimulus> provides the apertures as a cell vector of R x C x time.
%   values should be in [0,1].  the number of time points can differ across runs.
% <data> provides the data as a cell vector of voxels x time.  can also be
%   X x Y x Z x time.  the number of time points should match the number of 
%   time points in <stimulus>.
% <tr> is the TR in seconds (e.g. 1.5)
% <options> (optional) is a struct with the following fields:
%   <gaussianmode> (optional) is
%       'DoG' uses Difference of Gaussians model
%       'OG'  uses One Gaussian model
%       'GS'  uses Global Suppression model (One Gaussian with Global Suppression) 
%     default: 'DoG'.
%   <prfmodel> (optional) is
%       'CSS'     uses Compressive Spatial Summation model
%       'fixexpt' uses Compressive Spatial Summation model fixing <expt> at <typicalexpt> 
%       'linear'  uses linear model (same as <prfmodel='fixexpt'> with <typicalexpt=1.0>)
%     default: 'CSS'.
%   <vxs> (optional) is a vector of voxel indices to analyze.  (we automatically
%     sort the voxel indices and ensure uniqueness.)  default is 1:V where 
%     V is total number of voxels.  to reduce computational time, you may want 
%     to create a binary brain mask, perform find() on it, and use the result as <vxs>.
%   <wantglmdenoise> (optional) is whether to use GLMdenoise to determine
%     nuisance regressors to add into the PRF model.  note that in order to use
%     this feature, there must be at least two runs (and conditions must repeat
%     across runs).  we automatically determine the GLM design matrix based on
%     the contents of <stimulus>.  special case is to pass in the noise regressors 
%     directly (e.g. from a previous call).  default: 0.
%   <hrf> (optional) is a column vector with the hemodynamic response function (HRF)
%     to use in the model.  the first value of <hrf> should be coincident with the onset
%     of the stimulus, and the HRF should indicate the timecourse of the response to
%     a stimulus that lasts for one TR.  default is to use a canonical HRF (calculated
%     using getcanonicalhrf(tr,tr)').
%   <maxpolydeg> (optional) is a non-negative integer indicating the maximum polynomial
%     degree to use for drift terms.  can be a vector whose length matches the number
%     of runs in <data>.  default is to use round(L/2) where L is the number of minutes
%     in the duration of a given run.
%   <seedmode> (optional) is a vector consisting of one or more of the
%     following values (we automatically sort and ensure uniqueness):
%       0 means use generic large PRF seed
%       1 means use generic small PRF seed
%       2 means use best seed based on super-grid
%       3 means use pre-defined PRF seed (need options.initialseed)
%     default: [0 1 2].  a special case is to pass <seedmode> as -2.  this causes the
%     best seed based on the super-grid to be returned as the final estimate, thereby
%     bypassing the computationally expensive optimization procedure.  further notes
%     on this case are given below.
%   <dogseed> (optional) is a vector consisting of one or more of the
%     following values (we automatically sort and ensure uniqueness):
%       0 means use generic DoG seed or best DoG seed based on super-grid for seedmode=2 (Difference of Gaussians)
%       1 means use OG seed (One Gaussians)
%       2 means use generic SS seed (Slightly Suppression for DoG optimization)
%       3 means use generic GS seed or best GS seed based on super-grid for seedmode=2 (Global Suppression)
%     default: [0 1 2] for DoG model, [1] for OG model, [3] for GS model. 
%   <initialseed> (optional) is a X-params-channels matrix consisting of initial seeds. (only valid when seedmode=3)
%   <xvalmode> (optional) is
%     0 means just fit all the data
%     1 means N-fold cross-validation (first half of runs; second half of runs)
%     2 means N-fold cross-validation (first half of each run; second half of each run)
%     3 means N-fold cross-validation (first alternative runs; second alternative runs)
%     4 means N-fold cross-validation (first alternative of each run; second alternative of each run)
%     default: 0.  (note that we round when halving.)
%   <xvalfold> (optional) is number of fold in cross-validation.
%     N means N-fold cross-validation
%     1 means leave-one-outcross-validation
%     default: 2.  (two-fold cross-validation)
%   <numperjob> (optional) is
%     [] means to run locally (not on the cluster)
%     N where N is a positive integer indicating the number of voxels to analyze in each 
%       cluster job.  this option requires a customized computational setup!
%     default: [].
%   <maxiter> (optional) is the maximum number of iterations.  default: 500.
%   <display> (optional) is 'iter' | 'final' | 'off'.  default: 'iter'.
%   <typicalgain> (optional) is a typical value for the gain in each time-series.
%     default: 10.
%   <typicalexpt> (optional) is a typical value for the exponent in each time-series.
%     This option is ignored for  <prfmodel='linear'>.
%     default: 0.5 for <prfmodel='CSS'>, 0.05 for <prfmodel='fixexpt'>.
%   <forcebounds> (optional) is flag to force to consider bounds in optimization:
%     0 means ignore bounds
%     1 means set bounds in pRF eccentricity (refer boundsecc)
%     2 means set bounds in pRF eccentricity and exponent ([0,1])
%     default: 0.
%   <boundsecc> (optional) is how many times bounds of eccentricity is larger than stimulus:
%     if 2 (default), then the pRF center of each x and y axis is estimated in double area of stimulus
%   <noneggain> (optional) is flag to restrict gain to be non-negative:
%     default: true
%
% Analyze pRF data and return the results.
%
% The results structure contains the following fields:
% <ang> contains pRF angle estimates.  Values range between 0 and 360 degrees.
%   0 corresponds to the right horizontal meridian, 90 corresponds to the upper vertical
%   meridian, and so on.
% <ecc> contains pRF eccentricity estimates.  Values are in pixel units with a lower
%   bound of 0 pixels.
% <rfsize> contains pRF size estimates.  pRF size is defined as sigma/sqrt(n) where
%   sigma is the standard of the 2D Gaussian and n is the exponent of the power-law
%   function.  Values are in pixel units with a lower bound of 0 pixels.
% <expt> contains pRF exponent estimates.
% <gain> contains pRF gain estimates.  Values are in the same units of the data
%   and are constrained to be non-negative.
% <R2> contains R^2 values that indicate the goodness-of-fit of the model to the data.
%   Values are in percentages and generally range between 0% and 100%.  The R^2 values
%   are computed after projecting out polynomials from both the data and the model fit.
%   (Because of this projection, R^2 values can sometimes drop below 0%.)  Note that
%   if cross-validation is used (see <xvalmode>), the interpretation of <R2> changes
%   accordingly.
% <resnorms> and <numiters> contain optimization details (residual norms and 
%   number of iterations, respectively).
% <meanvol> contains the mean volume, that is, the mean of each voxel's time-series.
% <noisereg> contains a record of the noise regressors used in the model.
% <params> contains a record of the raw parameter estimates that are obtained internally
%   in the code.  These raw parameters are transformed to a more palatable format for
%   the user (as described above).
% <options> contains a record of the options used in the call to analyzePRF.m.
%
% Details on the pRF model:
% - Before analysis, we zero out any voxel that has a non-finite value or has all zeros
%   in at least one of the runs.  This prevents weird issues due to missing or bad data.
% - The pRF model that is fit is similar to that described in Dumoulin and Wandell (2008),
%   except that a static power-law nonlinearity is added to the model.  This new model, 
%   called the Compressive Spatial Summation (CSS) model, is described in Kay, Winawer, 
%   Mezer, & Wandell (2013).
% - The model involves computing the dot-product between the stimulus and a 2D isotropic
%   Gaussian, raising the result to an exponent, scaling the result by a gain factor,
%   and then convolving the result with a hemodynamic response function (HRF).  Polynomial
%   terms are included (on a run-by-run basis) to model the baseline signal level.
% - The 2D isotropic Gaussian is scaled such that the summation of the values in the
%   Gaussian is equal to one.  This eases the interpretation of the gain of the model.
% - The exponent parameter in the model is constrained to be non-negative.
% - The gain factor in the model is constrained to be non-negative; this aids the 
%   interpretation of the model (e.g. helps avoid voxels with negative BOLD responses
%   to the stimuli).
% - The workhorse of the analysis is fitnonlinearmodel.m, which is essentially a wrapper 
%   around routines in the MATLAB Optimization Toolbox.  We use the Levenberg-Marquardt 
%   algorithm for optimization, minimizing squared error between the model and the data.
% - A two-stage optimization strategy is used whereby all parameters excluding the
%   exponent parameter are first optimized (holding the exponent parameter fixed) and 
%   then all parameters are optimized (including the exponent parameter).  This 
%   strategy helps avoid local minima.
%
% Regarding GLMdenoise:
% - If the <wantglmdenoise> option is specified, we derive noise regressors using
%   GLMdenoise prior to model fitting.  This is done by creating a GLM design matrix
%   based on the contents of <stimulus> and then using this design matrix in conjunction
%   with GLMdenoise to analyze the data.  The noise regressors identified by GLMdenoise
%   are then used in the fitting of the pRF models (the regressors enter the model
%   additively, just like the polynomial regressors).
%
% Regarding seeding issues:
% - To minimize the impact of local minima, the default strategy is to perform full 
%   optimizations starting from three different initial seeds.
% - The first seed is a generic large pRF that is centered with respect to the stimulus,
%   has a pRF size equal to 1/4th of the stimulus extent (thus, +/- 2 pRF sizes matches
%   the stimulus extent), and has an exponent of 0.5.
% - The second seed is a generic small pRF that is just like the first seed except has
%   a pRF size that is 10 times smaller.
% - The third seed is a "supergrid" seed that is identified by performing a quick grid
%   search prior to optimization (similar in spirit to methods described in Dumoulin and 
%   Wandell, 2008).  In this procedure, a list of potential seeds is constructed by 
%   exploring a range of eccentricities, angles, and exponents.  For each potential 
%   seed, the model prediction is computed, and the seed that produces the closest 
%   match to the data is identified.  Note that the supergrid seed may be different
%   for different voxels.
%
% Regarding the "quick" mode:
% - When <seedmode> is -2, optimization is not performed and instead the best seed
%   based on the super-grid is returned as the final estimate.  If this case is used,
%   we automatically enforce that:
%   - opt.xvalmode is 0
%   - opt.vxs is []
%   - opt.numperjob is []
%   Also, in terms of outputs:
%   - The <gain> output is not estimated, and gain values are just returned as <typicalgain>.
%   - The <R2> output will contain correlation values (r) that range between -1 and 1.
%     These correlation values reflect the correlation between the model prediction and the
%     data after projecting out polynomial regressors and the noise regressors (if
%     <wantglmdenoise> is specified) from both the model prediction and the data.
%   - The <resnorms> and <numiters> outputs will be empty.
%
% history:
% 2023/06/28 - Yuasa: implement <seedmode> 3 to pass-through initial seeds.
% 2023/06/23 - Yuasa: implement N-fold cross-validation mode.
% 2022/12/05 - Yuasa: implement alternative cross-validation mode.
% 2022/05/16 - Yuasa: enable to modify bounds on eccentricity.
% 2020/02/25 - Yuasa: make bounds on eccentricity more strict.
%                     (triple of visual field -> double of visual field)
%                     separate <prfmodel='fixexpt'> from <prfmodel='linear'>
% 2020/02/02 - Yuasa: optimize gain seed in super-grid;
%                     super-grid search for GS seed.
% 2020/01/29 - Yuasa: enable to modify exponent seed by <typicalexpt>.
% 2020/01/27 - Yuasa: enable bounds with 'levenberg-marquardt' algorithm;
%                     add option to select OG/DoG and linear/CSS models.
% 2020/01/21 - Yuasa: output cross-validation accuracy.
% 2019/12/24 - Yuasa: add <dogseed> option.
% 2019/12/12 - Yuasa: increase initial seeds to avoid lower model fitting than OG model.
%                     this change requires much more computational time.
% 2019/12/06 - Yuasa: upate DoG model function.
% 2019/11/25 - Yuasa: change to estimate DOG model parameters.
% 2019/11/22 - Yuasa: modify model function for DoG model with fixed parameters.
% 2019/08/23 - Major change: the <stimulus> variable is now no longer forced to become
%              single format. This means the user controls whether computations are done
%              in double or single format. Please note that behavior (including finicky
%              local minima issues) can be highly dependent on the format. Parameter estimates
%              may be substantially more accurate (and may take substantially more computational
%              time / iterations to converge) if computations are performed in double format.
% 2015/02/07 - version 1.2
% 2015/02/07 - make analyzePRFcomputesupergridseeds.m less memory intensive
% 2014/11/10 - implement <wantglmdenoise> for the <seedmode> -2 case.
%              also, now the super-grid seed computation now
%              takes into account noise regressors (previously, the
%              noise regressors were ignored).
% 2014/10/20 - add -2 case for <seedmode>
% 2014/06/17 - version 1.1
% 2014/06/15 - add inputs <hrf> and <maxpolydeg>.
% 2014/06/10 - version 1.0
% 2014/04/27 - gain seed is now set to 0; add gain to the output
% 2014/04/29 - use typicalgain now (default 10). allow display input.

% internal notes:
% - for cluster mode, need to make sure fitnonlinearmodel is compiled (compilemcc.m)
% - to check whether local minima are a problem, can look at results.resnorms

% Dependency: setboundsinfunc, replelement

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% REPORT

fprintf('*** analyzePRF: started at %s. ***\n',datestr(now));
stime = clock;  % start time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INTERNAL CONSTANTS

% define
remotedir = '/scratch/knk/input/';
remotedir2 = '/scratch/knk/output/';
remotelogin = 'knk@login2.chpc.wustl.edu';
remoteuser = 'knk';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SETUP AND PREPARATION

% massage cell inputs
if ~iscell(stimulus)
  stimulus = {stimulus};
end
if ~iscell(data)
  data = {data};
end

% calc
is3d = size(data{1},4) > 1;
if is3d
  dimdata = 3;
  dimtime = 4;
  xyzsize = sizefull(data{1},3);
else
  dimdata = 1;
  dimtime = 2;
  xyzsize = size(data{1},1);
end
numvxs = prod(xyzsize);

% calc
res = sizefull(stimulus{1},2);
resmx = max(res);
numruns = length(data);

% deal with inputs
if ~exist('options','var') || isempty(options)
  options = struct();
end
if ~isfield(options,'gaussianmode') || isempty(options.gaussianmode)
  options.gaussianmode = 'dog';
else
  assert(ischar(options.gaussianmode),'gaussianmode is invalid');
  assert(ismember(lower(options.gaussianmode),{'dog','og','gs'}),'gaussianmode ''%s'' is invalid',options.gaussianmode);
  options.gaussianmode = lower(options.gaussianmode);
end
if ~isfield(options,'prfmodel') || isempty(options.prfmodel)
  options.prfmodel = 'css';
else
  assert(ischar(options.prfmodel),'prfmodel is invalid');
  assert(ismember(lower(options.prfmodel),{'css','fixexpt','linear'}),'prfmodel ''%s'' is invalid',options.prfmodel);
  options.prfmodel = lower(options.prfmodel);
end
if ~isfield(options,'vxs') || isempty(options.vxs)
  options.vxs = 1:numvxs;
end
if ~isfield(options,'wantglmdenoise') || isempty(options.wantglmdenoise)
  options.wantglmdenoise = 0;
end
if ~isfield(options,'hrf') || isempty(options.hrf)
  options.hrf = [];
end
if ~isfield(options,'maxpolydeg') || isempty(options.maxpolydeg)
  options.maxpolydeg = [];
end
if ~isfield(options,'seedmode') || isempty(options.seedmode)
  options.seedmode = [0 1 2];
end
if ~isfield(options,'dogseed') || isempty(options.dogseed)
  switch options.gaussianmode
    case 'dog', options.dogseed = [0 1 2];
    case 'og',  options.dogseed = 1;
    case 'gs',  options.dogseed = 3;
  end
end
if ~isfield(options,'initialseed') || isempty(options.initialseed)
    if all(options.seedmode == 3)
        error('''initialseed'' is required under seedmode = 3.');
    elseif ismember(3,options.seedmode)
        warning('''initialseed'' is not specified. Skip seedmode = 3.');
        options.seedmode(ismember(options.seedmode,3)) = [];
    end
end
if ~isfield(options,'xvalmode') || isempty(options.xvalmode)
  options.xvalmode = 0;
end
if ~isfield(options,'xvalfold') || isempty(options.xvalfold)
  options.xvalfold = 2;
end
if ~isfield(options,'numperjob') || isempty(options.numperjob)
  options.numperjob = [];
end
if ~isfield(options,'maxiter') || isempty(options.maxiter)
  options.maxiter = 500;
end
if ~isfield(options,'display') || isempty(options.display)
  options.display = 'iter';
end
if ~isfield(options,'typicalgain') || isempty(options.typicalgain)
  options.typicalgain = 10;
end
if ~isfield(options,'typicalexpt') || isempty(options.typicalexpt)
  switch options.prfmodel
    case 'css',     options.typicalexpt = 0.5;
    case 'fixexpt', options.typicalexpt = 0.05;
  end
end
switch options.prfmodel
    case 'linear',  options.typicalexpt = 1;
end
if ~isfield(options,'forcebounds') || isempty(options.forcebounds)
  options.forcebounds = 0;
end
if ~isfield(options,'boundsecc') || isempty(options.boundsecc)
  options.boundsecc = [-2 2];
elseif numel(options.boundsecc)==1
  options.boundsecc = [-1 1].*options.boundsecc;
end
if ~isfield(options,'noneggain') || isempty(options.noneggain)
  options.noneggain = true;
end

% massage
wantquick = isequal(options.seedmode,-2);
options.seedmode = union(options.seedmode(:),[]);

% massage more
if wantquick
  options.xvalmode = 0;
  options.vxs = 1:numvxs;
  options.numperjob = [];
  options.dogseed = 0;  % super-grid DoG seed
end

% calc
usecluster = ~isempty(options.numperjob);

% prepare stimuli
for p=1:length(stimulus)
  stimulus{p} = squish(stimulus{p},2)';  % frames x pixels
  stimulus{p} = [stimulus{p} p*ones(size(stimulus{p},1),1)];  % add a dummy column to indicate run breaks
%% REMOVED ON AUG 23 2019!  THIS CHANGES PAST BEHAVIOR.
%%  stimulus{p} = single(stimulus{p});  % make single to save memory
end

% deal with data badness (set bad voxels to be always all 0)
bad = cellfun(@(x) any(~isfinite(x),dimtime) | all(x==0,dimtime),data,'UniformOutput',0);  % if non-finite or all 0
bad = any(cat(dimtime,bad{:}),dimtime);  % badness in ANY run
for p=1:numruns
  data{p}(repmat(bad,[ones(1,dimdata) size(data{p},dimtime)])) = 0;
end

% calc mean volume
meanvol = mean(catcell(dimtime,data),dimtime);

% what HRF should we use?
if isempty(options.hrf)
  options.hrf = getcanonicalhrf(tr,tr)';
end
numinhrf = length(options.hrf);

% what polynomials should we use?
if isempty(options.maxpolydeg)
  options.maxpolydeg = cellfun(@(x) round(size(x,dimtime)*tr/60/2),data);
end
if isscalar(options.maxpolydeg)
  options.maxpolydeg = repmat(options.maxpolydeg,[1 numruns]);
end
fprintf('using the following maximum polynomial degrees: %s\n',mat2str(options.maxpolydeg));

% initialize cluster stuff
if usecluster
  localfilestodelete = {};
  remotefilestodelete = {};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURE OUT NOISE REGRESSORS

if isequal(options.wantglmdenoise,1)
  noisereg = analyzePRFcomputeGLMdenoiseregressors(stimulus,data,tr);
elseif isequal(options.wantglmdenoise,0)
  noisereg = [];
else
  noisereg = options.wantglmdenoise;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREPARE MODEL
eccrate = (options.boundsecc+1)./2;
lb = [1+eccrate(1).*res(1) 1+eccrate(1).*res(2) 0   0   0   1   0];
ub = [eccrate(2).*res(1)   eccrate(2).*res(2)   Inf Inf Inf Inf 1];
if options.forcebounds == 2 % upper limit of exponent
  ub(5) = 1;
end
if ~options.noneggain % lower limit of gain (also need to set in model)
  lb(4) = -inf;
end

% pre-compute some cache
[d,xx,yy] = makegaussian2d(resmx,2,2,2,2);

% define the model (parameters are R C S G N + SR GR) (embedded bounds: G>=0, N>=0)
modelfunbase = @(pp,dd) conv2run(modeldogcss(pp(1:5),pp(6:end),dd,res,xx,yy,0,0,options.noneggain),options.hrf,dd(:,prod(res)+1));
if ~options.forcebounds,	modelfun = @(pp,dd) modelfunbase([pp nan nan true],dd);
else,      modelfun = @(pp,dd)setboundsinfunc(lb,ub,[],[],@(pp,dd) modelfunbase([pp nan nan true],dd),pp,dd,'OutputMode','nan');
end
switch options.prfmodel
    case 'css'
        switch options.gaussianmode
            case 'dog'
model = {{[] [replelement(lb,5:7,NaN);
              ub                      ] modelfun} ...
         {@(ss)ss [replelement(lb,6:7,NaN);
                   ub                      ] @(ss)modelfun} ...
         {@(ss)ss [lb;
                   ub ] @(ss)modelfun}};
            case 'og'
model = {{[] [replelement(lb,5:7,NaN);
              ub                      ] modelfun} ...
         {@(ss)ss [replelement(lb,6:7,NaN);
                   ub                      ] @(ss)modelfun}};
            case 'gs'
model = {{[] [replelement(lb,5:7,NaN);
              ub                      ] modelfun} ...
         {@(ss)ss [replelement(lb,6:7,NaN);
                   ub                      ] @(ss)modelfun} ...
         {@(ss)ss [replelement(lb,6,NaN);
                   ub                    ] @(ss)modelfun}};
        end
    case {'fixexpt','linear'}
        switch options.gaussianmode
            case 'dog'
model = {{[] [replelement(lb,5:7,NaN);
              ub                      ] modelfun} ...
         {@(ss)ss [replelement(lb,5,NaN);
                   ub                  ] @(ss)modelfun}};
            case 'og'
model = {{[] [replelement(lb,5:7,NaN);
              ub                      ] modelfun}};
            case 'gs'
model = {{[] [replelement(lb,5:7,NaN);
              ub                      ] modelfun} ...
         {@(ss)ss [replelement(lb,5:6,NaN);
                   ub                      ] @(ss)modelfun}};
        end
end
               

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREPARE SEEDS

% init
seeds = [];
seedsOG  = [1 0];
seedsDoG = [2 0.5];
seedsSS  = [100 1e-3];
seedsGS  = [Inf 0.1];

% generic large seed
if ismember(0,options.seedmode)
  genseed = [(1+res(1))/2 (1+res(2))/2 resmx/4*sqrt(options.typicalexpt) options.typicalgain options.typicalexpt];
  if ismember(0,options.dogseed)
  seeds = [seeds;
           genseed seedsDoG];
  end
  if ismember(1,options.dogseed)
  seeds = [seeds;
           genseed seedsOG];
  end
  if ismember(2,options.dogseed)
  seeds = [seeds;
           genseed seedsSS];
  end
  if ismember(3,options.dogseed)
  seeds = [seeds;
           genseed seedsGS];
  end
end

% generic small seed
if ismember(1,options.seedmode)
  genseed = [(1+res(1))/2 (1+res(2))/2 resmx/4*sqrt(options.typicalexpt)/10 options.typicalgain options.typicalexpt];
  if ismember(0,options.dogseed)
  seeds = [seeds;
           genseed seedsDoG];
  end
  if ismember(1,options.dogseed)
  seeds = [seeds;
           genseed seedsOG];
  end
  if ismember(2,options.dogseed)
  seeds = [seeds;
           genseed seedsSS];
  end
  if ismember(3,options.dogseed)
  seeds = [seeds;
           genseed seedsGS];
  end
end

% super-grid seed
if any(ismember([2 -2],options.seedmode))
  datclass = class(data{1});
  data = cellfun(@single,data,'UniformOutput',false);   % temporally convert data to single
  if ismember(options.prfmodel,{'linear','fixexpt'}),  modelfun0 = @(pp,dd) modelfunbase([replelement(pp,5,options.typicalexpt) seedsOG],dd);
  else,                                                modelfun0 = @(pp,dd) modelfunbase([pp seedsOG],dd);
  end
  [ogsupergridseeds,ogrvalues] = analyzePRFcomputesupergridseeds(res,stimulus,data,modelfun0, ...
                                                   options.maxpolydeg,dimdata,dimtime, ...
                                                   options.typicalgain,noisereg);
  if ismember(options.prfmodel,{'linear','fixexpt'})  
      if is3d, ogsupergridseeds(:,:,:,5) = options.typicalexpt; 
      else,    ogsupergridseeds(:,5)     = options.typicalexpt; 
      end
  end
  [ogsupergridseeds]         = analyzePRFcomputegainseed(ogsupergridseeds,stimulus,data,modelfun0, ...
                                                   options.maxpolydeg,dimdata,dimtime, noisereg);
  if ismember(0,options.dogseed)
  [supergridseeds,rvalues] = analyzePRFdogcomputesupergridseeds(ogsupergridseeds,stimulus,data,modelfun, ...
                                                   options.maxpolydeg,dimdata,dimtime, ...
                                                   [],noisereg);
  else
   supergridseeds = [];
   rvalues        = [];
  end
  if ismember(1,options.dogseed)
   fullsupergridseeds        = cat(dimdata+1,ogsupergridseeds, ...
                               reshape(repmat(seedsOG,prod(sizefull(ogsupergridseeds,dimdata)),1),[sizefull(ogsupergridseeds,dimdata),2]));
   supergridseeds = cat(dimdata+2, supergridseeds, fullsupergridseeds);
   rvalues        = cat(dimdata+1, rvalues, ogrvalues);
  end
  if ismember(2,options.dogseed)
   fullsupergridseeds        = cat(dimdata+1,ogsupergridseeds, ...
                               reshape(repmat(seedsSS,prod(sizefull(ogsupergridseeds,dimdata)),1),[sizefull(ogsupergridseeds,dimdata),2]));
   supergridseeds = cat(dimdata+2, supergridseeds, fullsupergridseeds);
   rvalues        = cat(dimdata+1, rvalues, ogrvalues);
  end
  if ismember(3,options.dogseed)
  [fullsupergridseeds,fullrvalues] = analyzePRFgscomputesupergridseeds(ogsupergridseeds,stimulus,data,modelfun, ...
                                                   options.maxpolydeg,dimdata,dimtime, ...
                                                   [],noisereg);
   supergridseeds = cat(dimdata+2, supergridseeds, fullsupergridseeds);
   rvalues        = cat(dimdata+1, rvalues, fullrvalues);
  end
  data = cellfun(@(d) eval([datclass '(d)']),data,'UniformOutput',false);  % reconvert data
end

% pre-definedd seed
if ismember(3,options.seedmode)
  if ~exist('supergridseeds','var')
   supergridseeds = [];
  end
  if ndims(options.initialseed)<=3 && size(options.initialseed,3) == numvxs
      initialseed = permute(reshape(options.initialseed,[sizefull(options.initialseed,2),xyzsize]),[3:(2+dimdata),2,1]);
  elseif ndims(options.initialseed)<=3 && size(options.initialseed,1) == numvxs
      initialseed = reshape(options.initialseed,[xyzsize,size(options.initialseed,2:ndims(options.initialseed))]);
  elseif ndims(options.initialseed)>=dimdata && all(sizefull(options.initialseed,dimdata)==xyzsize)
      initialseed = options.initialseed;
  elseif ndims(options.initialseed)>=dimdata && all(size(options.initialseed,3:(2+dimdata))==xyzsize)
      initialseed = permute(options.initialseed,[3:(2+dimdata),2,1]);
  else
      error('The size of initialseed is invalid.');
  end
  if size(initialseed,dimdata+1)==5
    for p=1:size(initialseed,dimdata+2)
      if is3d, ogsupergridseeds  = initialseed(:,:,:,:,p);
      else,    ogsupergridseeds  = initialseed(:,:,p);
      end
      if ismember(0,options.dogseed)
      [fullsupergridseeds] = analyzePRFdogcomputesupergridseeds(ogsupergridseeds,stimulus,data,modelfun, ...
                                                       options.maxpolydeg,dimdata,dimtime, ...
                                                       [],noisereg);
       supergridseeds = cat(dimdata+2, supergridseeds, fullsupergridseeds);
      end
      if ismember(1,options.dogseed)
       fullsupergridseeds  = cat(dimdata+1,ogsupergridseeds, ...
                             reshape(repmat(seedsOG,prod(sizefull(ogsupergridseeds,dimdata)),1),[sizefull(ogsupergridseeds,dimdata),2]));
       supergridseeds = cat(dimdata+2, supergridseeds, fullsupergridseeds);
      end
      if ismember(2,options.dogseed)
       fullsupergridseeds  = cat(dimdata+1,ogsupergridseeds, ...
                             reshape(repmat(seedsSS,prod(sizefull(ogsupergridseeds,dimdata)),1),[sizefull(ogsupergridseeds,dimdata),2]));
       supergridseeds = cat(dimdata+2, supergridseeds, fullsupergridseeds);
      end
      if ismember(3,options.dogseed)
      [fullsupergridseeds] = analyzePRFgscomputesupergridseeds(ogsupergridseeds,stimulus,data,modelfun, ...
                                                       options.maxpolydeg,dimdata,dimtime, ...
                                                       [],noisereg);
       supergridseeds = cat(dimdata+2, supergridseeds, fullsupergridseeds);
      end
    end
  else
       supergridseeds = cat(dimdata+2, supergridseeds, initialseed);
  end

end

% make a function that individualizes the seeds
if exist('supergridseeds','var')
  seedfun = @(vx) [[seeds];
                   [permute(subscript(squish(supergridseeds,dimdata),{vx ':' ':'}),[3,2,1])]];
else
  seedfun = @(vx) [seeds];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PERFORM OPTIMIZATION

% if this is true, we can bypass all of the optimization stuff!
if wantquick

else

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREPARE RESAMPLING STUFF

  % define wantresampleruns and resampling
  switch options.xvalmode
  case 0
    wantresampleruns = [];
    resampling = 0;
  case 1
    wantresampleruns = 1;
    if options.xvalfold == 1,    xvalfold = length(data);
    else,                        xvalfold = options.xvalfold;
    end
    resamplingindex = floor(((1:length(data))-1)./length(data).*xvalfold)+1;
    resampling = ones(max(resamplingindex),length(data));
    for nx=1:max(resamplingindex)
        resampling(nx,resamplingindex==(max(resamplingindex)-nx+1)) = -1;
    end
  case 2
    wantresampleruns = 0;
    resampling = [];
    for p=1:length(data)
      if options.xvalfold == 1,    xvalfold = size(data{p},2);
      else,                        xvalfold = options.xvalfold;
      end
      resamplingindex = floor(((1:size(data{p},2))-1)./size(data{p},2).*xvalfold)+1;
      resamplingP = ones(max(resamplingindex),size(data{p},2));
      for nx=1:max(resamplingindex)
          resamplingP(nx,resamplingindex==(max(resamplingindex)-nx+1)) = -1;
      end
      resampling = cat(2,resampling,resamplingP);
    end
  case 3
    wantresampleruns = 1;
    if options.xvalfold == 1,    xvalfold = length(data);
    else,                        xvalfold = options.xvalfold;
    end 
    resamplingindex = mod((1:length(data))-1,xvalfold)+1;
    resampling = ones(max(resamplingindex),length(data));
    for nx=1:max(resamplingindex)
        resampling(nx,resamplingindex==(max(resamplingindex)-nx+1)) = -1;
    end
  case 4
    wantresampleruns = 0;
    resampling = [];
    for p=1:length(data)
      if options.xvalfold == 1,    xvalfold = size(data{p},2);
      else,                        xvalfold = options.xvalfold;
      end
      resamplingindex = mod((1:size(data{p},2))-1,xvalfold)+1;
      resamplingP = ones(max(resamplingindex),size(data{p},2));
      for nx=1:max(resamplingindex)
          resamplingP(nx,resamplingindex==(max(resamplingindex)-nx+1)) = -1;
      end
      resampling = cat(2,resampling,resamplingP);
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREPARE STIMULUS AND DATA

  %%%%% CLUSTER CASE

  if usecluster

    % save stimulus and transport to the remote server
    while 1
      filename0 = sprintf('stim%s.mat',randomword(5));  % file name
      localfile0 = [tempdir '/' filename0];             % local path to file
      remotefile0 = [remotedir '/' filename0];          % remote path to file
  
      % redo if file already exists locally or remotely
      if exist(localfile0) || 0==unix(sprintf('ssh %s ls %s',remotelogin,remotefile0))
        continue;
      end
  
      % save file and transport it
      save(localfile0,'stimulus');
      assert(0==unix(sprintf('rsync -av %s %s:"%s/"',localfile0,remotelogin,remotedir)));
  
      % record
      localfilestodelete{end+1} = localfile0;
      remotefilestodelete{end+1} = remotefile0;
  
      % stop
      break;
    end
    clear stimulus;  % don't let it bleed through anywhere!

    % define stimulus
    stimulus = @() loadmulti(remotefile0,'stimulus');

    % save data and transport to the remote server
    while 1
      filename0 = sprintf('data%s',randomword(5));   % directory name that will contain 001.bin, etc.
      localfile0 = [tempdir '/' filename0];          % local path to dir
      remotefile0 = [remotedir '/' filename0];       % remote path to dir
  
      % redo if dir already exists locally or remotely
      if exist(localfile0) || 0==unix(sprintf('ssh %s ls %s',remotelogin,remotefile0))
        continue;
      end
  
      % save files and transport them
      assert(mkdir(localfile0));
      for p=1:numruns
        savebinary([localfile0 sprintf('/%03d.bin',p)],'single',squish(data{p},dimdata)');  % notice squish
      end
      assert(0==unix(sprintf('rsync -av %s %s:"%s/"',localfile0,remotelogin,remotedir)));

      % record
      localfilestodelete{end+1} = localfile0;
      remotefilestodelete{end+1} = remotefile0;

      % stop
      break;
    end
    clear data;

    % define data
    binfiles = cellfun(@(x) [remotefile0 sprintf('/%03d.bin',x)],num2cell(1:numruns),'UniformOutput',0);
    data = @(vxs) cellfun(@(x) double(loadbinary(x,'single',[0 numvxs],-vxs)),binfiles,'UniformOutput',0);

    % prepare the output directory name
    while 1
      filename0 = sprintf('prfresults%s',randomword(5));
      localfile0 = [tempdir '/' filename0];
      remotefile0 = [remotedir2 '/' filename0];
      if exist(localfile0) || 0==unix(sprintf('ssh %s ls %s',remotelogin,remotefile0))
        continue;
      end
      localfilestodelete{end+1} = localfile0;
      localfilestodelete{end+1} = [localfile0 '.mat'];  % after consolidation
      remotefilestodelete{end+1} = remotefile0;
      break;
    end
    outputdirlocal = localfile0;
    outputdirremote = remotefile0;
    outputdir = outputdirremote;

  %%%%% NON-CLUSTER CASE

  else

    stimulus = {stimulus};
    data = @(vxs) cellfun(@(x) subscript(squish(x,dimdata),{vxs ':'})',data,'UniformOutput',0);
    outputdir = [];

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREPARE OPTIONS

  % last-minute prep
  if iscell(noisereg)
    noiseregINPUT = {noisereg};
  else
    noiseregINPUT = noisereg;
  end

  % construct the options struct    % Algorithm: 'levenberg-marquardt', 'trust-region-reflective'
  opt = struct( ...
    'outputdir',outputdir, ...
    'stimulus',stimulus, ...
    'data',data, ...
    'vxs',options.vxs, ...
    'model',{model}, ...
    'seed',seedfun, ...
    'optimoptions',{{'Display' options.display 'Algorithm' 'levenberg-marquardt' 'MaxIter' options.maxiter}}, ...
    'wantresampleruns',wantresampleruns, ...
    'resampling',resampling, ...
    'metric',@calccod, ...
    'maxpolydeg',options.maxpolydeg, ...
    'wantremovepoly',1, ...
    'extraregressors',noiseregINPUT, ...
    'wantremoveextra',0, ...
    'dontsave',{{'modelfit' 'opt' 'vxsfull' 'modelpred' 'testdata'}});  % 'resnorms' 'numiters' 

          %  'outputfcn',@(a,b,c,d) pause2(.1) | outputfcnsanitycheck(a,b,c,1e-6,10) | outputfcnplot(a,b,c,1,d), ...
          %'outputfcn',@(a,b,c,d) pause2(.1) | outputfcnsanitycheck(a,b,c,1e-6,10) | outputfcnplot(a,b,c,1,d));
          %   % debugging:
          %   chpcstimfile = '/stone/ext1/knk/HCPretinotopy/conimagesB.mat';
          %   chpcdatadir2 = outputdir2;  % go back
          %   opt.outputdir='~/temp1';
          %   profile on;
          %   results = fitnonlinearmodel(opt,100,100);
          %   results = fitnonlinearmodel(opt,1,715233);
          %   profsave(profile('info'),'~/inout/profile_results');
          % %   modelfit = feval(modelfun,results.params,feval(stimulusINPUT));
          % %   thedata = feval(dataINPUT,52948);
          % %   pmatrix = projectionmatrix(constructpolynomialmatrix(304,0:3));
          % %   figure; hold on;
          % %   plot(pmatrix*thedata,'k-');
          % %   plot(pmatrix*modelfit,'r-');

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIT MODEL

  %%%%% CLUSTER CASE

  if usecluster

    % submit jobs
    jobnames = {};
    jobnames = [jobnames {makedirid(opt.outputdir,1)}];
    jobids = [];
    jobids = [jobids chpcrun(jobnames{end},'fitnonlinearmodel',options.numperjob, ...
                             1,ceil(length(options.vxs)/options.numperjob),[], ...
                             {'data' 'stimulus' 'bad' 'd' 'xx' 'yy' 'modelfun' 'model'})];

    % record additional files to delete
    for p=1:length(jobnames)
      remotefilestodelete{end+1} = sprintf('~/sgeoutput/job_%s.*',jobnames{p});  % .o and .e files
      remotefilestodelete{end+1} = sprintf('~/mcc/job_%s.mat',jobnames{p});
      localfilestodelete{end+1} = sprintf('~/mcc/job_%s.mat',jobnames{p});
    end

    % wait for jobs to finish
    sgewaitjobs(jobnames,jobids,remotelogin,remoteuser);

    % download the results
    assert(0==unix(sprintf('rsync -a %s:"%s" "%s/"',remotelogin,outputdirremote,tempdir)));

    % consolidate the results
    fitnonlinearmodel_consolidate(outputdirlocal);

    % load the results
    a1 = load([outputdirlocal '.mat']);

  %%%%% NON-CLUSTER CASE

  else

    a1 = fitnonlinearmodel(opt);

  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREPARE OUTPUT

% depending on which analysis we did (quick or full optimization),
% we have to get the outputs in a common format
if wantquick
  [rA,maxrA] = max(squish(rvalues,dimdata)',[],1,'linear');   % fits x voxels
  paramsA = permute(squish(supergridseeds,dimdata),[2 3 1]);  % parameters x fits x voxels
  paramsA = permute(paramsA(:,maxrA),[3 2 1]);                % fits x parameters x voxels
else
  paramsA = a1.params;                                        % fits x parameters x voxels
  rA = a1.trainperformance;                                   % fits x voxels
end

% calc
numfits = size(paramsA,1);

% init
clear results;
results.ang =      NaN*zeros(numvxs,numfits);
results.ecc =      NaN*zeros(numvxs,numfits);
results.expt =     NaN*zeros(numvxs,numfits);
results.rfsize =   NaN*zeros(numvxs,numfits);
results.R2 =       NaN*zeros(numvxs,numfits);
results.gain =     NaN*zeros(numvxs,numfits);
results.xR2 =      NaN*zeros(numvxs,numfits);
results.xval =     NaN*zeros(numvxs,1);
results.resnorms = cell(numvxs,1);
results.numiters = cell(numvxs,1);

% massage model parameters for output and put in 'results' struct
results.ang(options.vxs,:) =    permute(mod(atan2((1+res(1))/2 - paramsA(:,1,:), ...
                                                  paramsA(:,2,:) - (1+res(2))/2),2*pi)/pi*180,[3 1 2]);
results.ecc(options.vxs,:) =    permute(sqrt(((1+res(1))/2 - paramsA(:,1,:)).^2 + ...
                                             (paramsA(:,2,:) - (1+res(2))/2).^2),[3 1 2]);
results.expt(options.vxs,:) =   permute(posrect(paramsA(:,5,:)),[3 1 2]);
results.rfsize(options.vxs,:) = permute(abs(paramsA(:,3,:)) ./ sqrt(posrect(paramsA(:,5,:))),[3 1 2]);
results.R2(options.vxs,:) =     permute(rA,[2 1]);
results.gain(options.vxs,:) =   permute(paramsA(:,4,:),[3 1 2]);
if options.noneggain
results.gain(options.vxs,:) =   posrect(results.gain(options.vxs,:));
end
if ~wantquick
  results.resnorms(options.vxs) = a1.resnorms;
  results.numiters(options.vxs) = a1.numiters;
  if options.xvalmode > 0
      results.xR2(options.vxs,:) =     permute(a1.testperformance,[2 1]);
      results.xval(options.vxs)  =     a1.aggregatedtestperformance;
  end
end

% reshape
results.ang =      reshape(results.ang,      [xyzsize numfits]);
results.ecc =      reshape(results.ecc,      [xyzsize numfits]);
results.expt =     reshape(results.expt,     [xyzsize numfits]);
results.rfsize =   reshape(results.rfsize,   [xyzsize numfits]);
results.R2 =       reshape(results.R2,       [xyzsize numfits]);
results.gain =     reshape(results.gain,     [xyzsize numfits]);
results.xR2 =      reshape(results.xR2,      [xyzsize numfits]);
results.xval =     reshape(results.xval,     [xyzsize 1]);
results.resnorms = reshape(results.resnorms, [xyzsize 1]);
results.numiters = reshape(results.numiters, [xyzsize 1]);

% add some more stuff
results.meanvol =  meanvol;
results.noisereg = noisereg;
results.params =   paramsA;
results.options = options;

% save 'results' to a temporary file so we don't lose these precious results!
file0 = [tempname '.mat'];
fprintf('saving results to %s (just in case).\n',file0);
save(file0,'results');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CLEAN UP

% no clean up necessary in the quick case
if ~wantquick

  %%%%% CLUSTER CASE

  if usecluster

    % delete local files and directories
    for p=1:length(localfilestodelete)
      if exist(localfilestodelete{p},'dir')  % first dir, then file
        rmdir(localfilestodelete{p},'s');
      elseif exist(localfilestodelete{p},'file')
        delete(localfilestodelete{p});
      end
    end

    % delete remote files and directories
    for p=1:length(remotefilestodelete)
      assert(0==unix(sprintf('ssh %s "rm -rf %s"',remotelogin,remotefilestodelete{p})));
    end

  %%%%% NON-CLUSTER CASE

  else

  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% REPORT

fprintf('*** analyzePRF: ended at %s (%.1f minutes). ***\n', ...
        datestr(now),etime(clock,stime)/60);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% JUNK

% % define the model (parameters are R C S G N [HRF])
% modelfun = @(pp,dd) conv2run(posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res),makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*abs(pp(3))^2))); 0]) .^ posrect(pp(5)),pp(5+(1:numinhrf))',dd(:,prod(res)+1));
% model = {{[] [1-res(1)+1 1-res(2)+1 0    0   NaN repmat(NaN,[1 numinhrf]);
%               2*res(1)-1 2*res(2)-1 Inf  Inf Inf repmat(Inf,[1 numinhrf])] modelfun} ...
%          {@(ss)ss [1-res(1)+1 1-res(2)+1 0    0   0   repmat(NaN,[1 numinhrf]);
%                    2*res(1)-1 2*res(2)-1 Inf  Inf Inf repmat(Inf,[1 numinhrf])] @(ss)modelfun} ...
%          {@(ss)ss [1-res(1)+1 1-res(2)+1 0    0   0   repmat(-Inf,[1 numinhrf]);
%                    2*res(1)-1 2*res(2)-1 Inf  Inf Inf repmat(Inf,[1 numinhrf])] @(ss)modelfun}};
% 
% % if not fitting the HRF, exclude the last model step
% if ~wantfithrf
%   model = model(1:2);
% end
%wantfithrf = 0;    % for now, leave at 0

% results.hrf =      NaN*zeros(numvxs,numinhrf,numfits);
% results.hrf(options.vxs,:,:) =  permute(a1.params(:,5+(1:numinhrf),:),[3 2 1]);
% results.hrf =      reshape(results.hrf,      [xyzsize numinhrf numfits]);
