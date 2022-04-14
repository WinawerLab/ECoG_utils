function [seeds,rvalues] = analyzePRFdogcomputesupergridseeds(seeds,stimulus,data,modelfun,maxpolydeg,dimdata,dimtime,typicalgain,noisereg)

% function [seeds,rvalues] = analyzePRFdogcomputesupergridseeds(seeds,stimulus,data,modelfun,maxpolydeg,dimdata,dimtime,typicalgain,noisereg)
%
% <seeds> is is outputs of analyzePRFcomputesupergridseeds
% <stimulus> is a cell vector of time x (pixels+1)
% <data> is a cell vector of X x Y x Z x time (or XYZ x time)
% <modelfun> is a function that accepts parameters (pp) and stimuli (dd) and outputs predicted time-series (time x 1)
% <maxpolydeg> is a vector of degrees (one element for each run)
% <dimdata> is number of dimensions that pertain to voxels
% <dimtime> is the dimension that is the time dimension
% <typicalgain> is [] or a typical value for the gain in each time-series
%               if [] then typicalgain is copied from <seeds>
% <noisereg> is [] or a set of noise regressors (cell vector of matrices)
%
% this is an internal function called by analyzePRFdog.m.  this function returns <seeds>
% as a matrix of dimensions X x Y x Z x parameters (or XYZ x parameters)
% with the best seed from the super-grid.  also, returns <rvalues> as X x Y x Z
% (or XYZ x 1) with the corresponding correlation (r) values.
%
% history:
% 2020/01/31 - Yuasa: enable to pass-through gain seeds
% 2019/12/12 - Yuasa: modify for DoG
% 2015/02/07 - make less memory intensive

% internal notes:
% - note that the gain seed is fake (it is not set the correct value but instead
%   to the <typicalgain>)

% internal constants
srs = [(10.^linspace(1,2,10))/10 25 50 100 1000];  srs(1) = [];
grs = (10.^linspace(-6,-0.1,14));
% remove sigma=1, gain=0
% (initiate with gain=0 increases the risk of local minimum problem into OG)

% calc
numvxs = prod(sizefull(data{1},dimdata));  % total number of voxels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% calculate a long list of potential seeds

% calc
% construct full list of seeds (seeds x params) [SR GR]
fprintf('constructing seeds.\n');
[X,Y] = meshgrid(srs,grs);
allseeds = reshape(cat(3,X,Y),[],2);                % not include simga=1, gain=0
% allseeds = cat(1,[1 0],reshape(cat(3,X,Y),[],2));   % add simga=1, gain=0
cnt = size(allseeds,1);
% construct full list of seeds (seeds x params) [R C S G N]
seedslist = unique(squish(seeds,dimdata),'rows','stable');
X = reshape(permute(repmat(seedslist,[1,1,cnt]),[3,1,2]),[],size(seedslist,2));
Y = repmat(allseeds,size(seedslist,1),1);
allseeds = cat(2,X,Y);
% set gain as 1
allseeds(:,4) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% generate the predicted time-series for each seed

% generate predicted time-series [note that some time-series are all 0]
predts = zeros(sum(cellfun(@(x) size(x,1),stimulus)),size(allseeds,1),'single');  % time x seeds
temp = catcell(1,stimulus);
fprintf('generating super-grid time-series...'); tic
parfor p=1:size(allseeds,1)
  predts(:,p) = modelfun(allseeds(p,:),temp);
end
fprintf('done.'); toc
clear temp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% prepare data and model predictions

% construct polynomials, noise regressors, and projection matrix
pregressors = {};
for p=1:length(maxpolydeg)
  pregressors{p} = constructpolynomialmatrix(size(data{p},dimtime),0:maxpolydeg(p));
  if ~isempty(noisereg)
    pregressors{p} = cat(2,pregressors{p},noisereg{p});
  end
end
pmatrix = projectionmatrix(blkdiag(pregressors{:}));

% project out and scale to unit length
predts = unitlength(pmatrix*predts,                                1,[],0);  % time x seeds   [NOTE: some are all NaN]
  % OLD: datats = unitlength(pmatrix*squish(catcell(dimtime,data),dimdata)',1,[],0);  % time x voxels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% find the seed with the max correlation

% compute correlation and find maximum for each voxel
chunks = chunking(1:numvxs,100);
rvalues = {};
bestseedix = {};
fprintf('finding best seed for each voxel.\n');
parfor p=1:length(chunks)

  % project out and scale to unit length
  datats = unitlength(pmatrix*catcell(2,cellfun(@(x) subscript(squish(x,dimdata),{chunks{p} ':'}),data,'UniformOutput',0))',1,[],0);  % time x voxels

  % voxels x 1 with index of the best seed (max corr)
  [rvalues{p},bestseedix{p}] = max(datats'*predts,[],2);  % voxels x seeds -> max corr along dim 2 [NaN is ok]

end
rvalues = catcell(1,rvalues);        % voxels x 1
bestseedix = catcell(1,bestseedix);  % voxels x 1

% prepare output
if isempty(typicalgain), typicalgain = seeds(:,4);  end
rvalues = reshape(rvalues,[sizefull(data{1},dimdata) 1]);
seeds = allseeds(bestseedix,:);  % voxels x parameters
seeds(:,4) = typicalgain;        % set gain to typical gain
seeds = reshape(seeds,[sizefull(data{1},dimdata) size(allseeds,2)]);



%  predts(:,p) = modelfun([allseeds(p,:) flatten(hrf)],temp);
% <hrf> is T x 1 with the HRF
