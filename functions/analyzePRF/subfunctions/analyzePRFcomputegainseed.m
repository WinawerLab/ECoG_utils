function [seeds,rvalues] = analyzePRFcomputegainseed(seeds,stimulus,data,modelfun,maxpolydeg,dimdata,dimtime,noisereg)

% function [seeds,rvalues] = analyzePRFcomputegainseed(seeds,stimulus,data,modelfun,maxpolydeg,dimdata,dimtime,noisereg)
%
% <seeds> is is outputs of analyzePRFcomputesupergridseeds
% <stimulus> is a cell vector of time x (pixels+1)
% <data> is a cell vector of X x Y x Z x time (or XYZ x time)
% <modelfun> is a function that accepts parameters (pp) and stimuli (dd) and outputs predicted time-series (time x 1)
% <maxpolydeg> is a vector of degrees (one element for each run)
% <dimdata> is number of dimensions that pertain to voxels
% <dimtime> is the dimension that is the time dimension
% <noisereg> is [] or a set of noise regressors (cell vector of matrices)
%
% this is an internal function called by analyzePRFdog.m.  this function
% updates <gain> in <seeds>.
%
% history:
% 2020/02/02 - Yuasa: estimate typical gain for each voxel
% 2015/02/07 - make less memory intensive

% internal notes:
% - note that the gain seed is fake (it is not set the correct value but instead
%   to the <typicalgain>)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% calculate a long list of potential seeds

% calc
allseeds = squish(seeds,dimdata);  % chop because of the omission above
% set gain as 1
allseeds(:,4) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% generate the predicted time-series for each seed

% generate predicted time-series [note that some time-series are all 0]
predts = zeros(sum(cellfun(@(x) size(x,1),stimulus)),size(allseeds,1),class(data{1}));  % time x seeds
temp = catcell(1,stimulus);
fprintf('generating time-series in each voxel...'); tic
parfor p=1:size(allseeds,1)
  predts(:,p) = modelfun(allseeds(p,:),temp);
end
fprintf('done.'); toc
clear temp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% prepare data and model predictions

% construct polynomials, noise regressors, and projection matrix
pregressors = {};
for p=1:length(maxpolydeg)
  if isnan(maxpolydeg(p))
    pregressors{p} = zeros(size(data{p},dimtime),0);
  else
    pregressors{p} = constructpolynomialmatrix(size(data{p},dimtime),0:maxpolydeg(p));
  end
  if ~isempty(noisereg)
    pregressors{p} = cat(2,pregressors{p},noisereg{p});
  end
end
pmatrix = projectionmatrix(blkdiag(pregressors{:}));

% project out and get scale
predts   = pmatrix*predts;  predts(predts<0) = 0;
[predts, predgain] = unitlength(predts, 1,[],0);  % time x seeds   [NOTE: some are all NaN]
  % OLD: datats = unitlength(pmatrix*squish(catcell(dimtime,data),dimdata)',1,[],0);  % time x voxels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% find the seed with the max correlation

% compute data time series and get scale
  % project out and scale to unit length
  datats            = pmatrix*catcell(2,cellfun(@(x) squish(x,dimdata),data,'UniformOutput',0))';
  datats(datats<0)  = 0;
  [datats, datagain] = unitlength(datats,1,[],0);  % time x voxels

% prepare output
gainseeds = datagain./predgain;
seeds(:,4) = gainseeds;        % set gain to typical gain
seeds = reshape(seeds,[sizefull(data{1},dimdata) size(allseeds,2)]);





%  predts(:,p) = modelfun([allseeds(p,:) flatten(hrf)],temp);
% <hrf> is T x 1 with the HRF
