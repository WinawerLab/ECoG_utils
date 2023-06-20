function [modelts, datats] = reconPRFdog(result,stimulus,data,opts)

% [modelts] = reconPRFdog(result,stimulus,option)
% [modelts] = reconPRFdog(params,stimulus,option)
%   returns model time series estimated from pRF parametes in result and stimulus.
% [modelts, datats] = reconPRFdog(result,stimulus,data,option)
% [modelts, datats] = reconPRFdog(result,datastcuct,option)
%   returns model and data time series estimated from pRF parametes in result and stimulus.
% 
% result:           output of analyzePRF
% params:           params field of the output of analyzePRF
% stimulus:         time series of the stimulus apertures as a cell array of R x C x time.
% data:             time series of the data as a cell array of voxels x time
% option:           structure with following fields.
%   whichelec:      [](default), cell array of channel names, # of channels, logical indices
%   stimres:        [](default), cell array of channel names, # of channels, logical indices
%   flipgain:       'yes', 'no'(default)  % apply -1 to the gain
%   skipprojection: 'yes', 'no'(default)  % apply polynomial projection
%   catrun:         'yes', 'no'(default)  % concatenate cell array across runs
%   catiter:        'yes', 'no'(default)  % concatenate cell array across iterations
%   catchan:        'yes', 'no'(default)  % concatenate cell array across channels
% 
%   Following options refer to result.options, but you can specify in option:
%       option.hrf
%       option.maxpolydeg
%       option.gaussianmode
% 
% modelts:          output time series of the model as a cell array of {chans}{iter x runs}(time x 1).
% datats:           output time series of the data as a cell array of {chans}{iter x runs}(time x 1).
%                   if catrun and catiter are 'yes':  {chans}(timeXruns x iter)
%                   if catchan is 'yes' besides:      (timeXruns x chans)

% Dependency: cellstrfind

% 20210723 Yuasa
% 20211122 Yuasa - update 'og' model
% 20221006 Yuasa - add options
% 20230620 Yuasa - delete unappropriate option

%%
%%% Check input
narginchk(2,4);
if nargin < 4,  opts = [];      end
if nargin < 3,  data = [];
elseif isstruct(data)&&~isfield(data,'datats')
                opts = data;
                data = [];
end
if isstruct(stimulus)
    if isempty(data)
        data     = stimulus;
        stimulus = [];
    else
        stimulus = stimulus.stimulus;
    end
end
if isempty(data), nargoutchk(0,1);
else,             nargoutchk(0,2);
end
if isempty(stimulus)
    if isfield(opts,'stimulus')&&~isempty(opts.stimulus)
        stimulus = opts.stimulus;
    elseif isfield(data,'stimulus')&&~isempty(data.stimulus)
        stimulus = data.stimulus;
    else
        error('stimulus is invalid');
    end
end
if isfield(opts,'channels')&&~isempty(opts.channels)
    channels = opts.channels;
elseif isfield(result,'channels')&&~isempty(result.channels)
    channels = result.channels;
elseif isfield(data,'channels')&&~isempty(data.channels)
    channels = data.channels;
else
    channels = '';
end
%-- correct data structure
if ~isstruct(result),   result = struct('params',result);   end
if isstruct(data),      data   = data.datats;               end

%%% Set options
SetDefault('result.options.hrf',1,0);
SetDefault('opts.hrf',result.options.hrf,0);
SetDefault('result.options.maxpolydeg',1,0);
SetDefault('opts.maxpolydeg',result.options.maxpolydeg,0);
SetDefault('result.options.gaussianmode','og',0);
SetDefault('opts.gaussianmode',result.options.gaussianmode,0);
SetDefault('opts.whichelec','*',0);
SetDefault('opts.flipgain','no',0);
SetDefault('opts.skipprojection','no',0);
SetDefault('opts.catrun','no',0);
SetDefault('opts.catiter','no',0);
SetDefault('opts.catchan','no',0);
SetDefault('opts.stimres',[],1);

%%% Set data & stimulus sizes (duplicate)
numvxs  = size(result.params,3);
numitr  = size(result.params,4);
if ~iscell(stimulus),      stimulus = {stimulus}; end
if ~iscell(stimulus{1}),   stimulus = {stimulus};  end
if length(stimulus)==1,    stimulus = repmat(stimulus,numvxs,1);    end
if isempty(data)
    [~,maxrunidx] = max(cellfun(@length,stimulus)); 
    data = cellfun(@(C) zeros(numvxs,size(C,3)),stimulus{maxrunidx},'UniformOutput',false);
    data = repmat(data,numitr,1);
else
    if ~iscell(data), data = {data}; end
end
numruns = size(data,2);

%%% Select channels
if islogical(opts.whichelec)
    whichelec = find(opts.whichelec);
elseif isnumeric(opts.whichelec)
    whichelec = opts.whichelec;
elseif strcmp(opts.whichelec,'*') && (isempty(channels)||~isfield(summary(channels),'name'))
    whichelec = 1:size(result.params,3);
else
    whichelec = cellstrfind(channels.name, opts.whichelec,0);
end
numelec = length(whichelec);

%%% Set parametes
hrf           = opts.hrf;
degs          = opts.maxpolydeg;
gaussianmode  = lower(opts.gaussianmode);
gainsign      = (-1)^intprtopt(opts.flipgain);
skipproj      = intprtopt(opts.skipprojection);
catrun        = intprtopt(opts.catrun) | numruns == 1;
catiter       = intprtopt(opts.catiter) | numitr == 1;
catchan       = intprtopt(opts.catchan) | numelec == 1;
res           = opts.stimres;
%-- projection
if skipproj, degs(:) = -1;  end
if numel(degs)<numruns
    assert(numel(degs)==1,'maxpolydeg option is invalid.');
    degs = repmat(degs,1,numruns);
end

%% Prepare model computation
%-- Get data length
ntime = zeros(numvxs,numruns);
for pp=1:numruns
    ntime(:,pp) = sum(~isnan(data{1,pp}),2);
end
%-- Set stimulus (change size)
if isempty(res)
    res = sizefull(stimulus{1}{1},2);
else
    if isscalar(res),       res = [res res];
    elseif ~isvector(res),  res = res(1:2);
    end
    stimulus = cellfun(...
            @(C) cellfun(@(S) imresize(S, res, 'nearest'),C,'UniformOutput',false),...
            stimulus,'UniformOutput',false);
end
resmx = max(res);

%%% Model computation
%-- Pre-compute cache for faster execution
[~,xx,yy] = makegaussian2d(resmx,2,2,2,2);

%-- Prepare the stimuli for use in the model
stimulusPP = repmat({{}},numvxs,1);
for vxs=1:numvxs
for pp=1:numruns
  stimulusPP{vxs}{pp} = squish(stimulus{vxs}{pp},2)';  % this flattens the image so that the dimensionality is now frames x pixels
  stimulusPP{vxs}{pp} = [stimulusPP{vxs}{pp} pp*ones(size(stimulusPP{vxs}{pp},1),1)];  % this adds a dummy column to indicate run breaks
end
end

%-- Set model
switch gaussianmode
    case {'dog','gs','lfs'}
        modelfun = @(pp,dd) conv2run(modeldogcss(pp(1:5),pp(6:end),dd,res,xx,yy,0,0),hrf,dd(:,prod(res)+1));
    case {'og'}
        modelfun = @(pp,dd) conv2run(posrect(pp(4)) * amppow((dd*[vflatten(placematrix(zeros(res),makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*abs(pp(3))^2))); 0]), posrect(pp(5))),hrf,dd(:,prod(res)+1));
end

%-- Construct projection matrices
polymatrix = repmat({{}},numvxs,1);
for vxs=1:numvxs
for pp=1:numruns
    polymatrix{vxs}{pp} = projectionmatrix(constructpolynomialmatrix(ntime(vxs,pp),0:degs(pp)));
end
end

%% Contsruct time series
datats  = repmat({{}},numvxs,1);
modelts = repmat({{}},numvxs,1);
for vxs=reshape(whichelec,1,[])
for pp=1:numruns
for itr=1:numitr
  nanplace    = isnan(data{itr,pp}(vxs,:));
  datats{vxs}{itr,pp}  = gainsign.*polymatrix{vxs}{pp}*data{itr,pp}(vxs,~nanplace)';
  modelts{vxs}{itr,pp} = gainsign.*polymatrix{vxs}{pp}*modelfun(result.params(1,:,vxs,itr),stimulusPP{vxs}{pp}(~nanplace,:));
end
end
end
datats  = datats(whichelec);
modelts = modelts(whichelec);

if catrun && catiter
    datats  = cellfun(@(C) cell2mat(permute(C,[2,3,1])),datats,'UniformOutput',false);
    modelts = cellfun(@(C) cell2mat(permute(C,[2,3,1])),modelts,'UniformOutput',false);
elseif catrun
    datats  = cellfun(@(C) cellcatcolumn(C,1),datats,'UniformOutput',false);
    modelts = cellfun(@(C) cellcatcolumn(C,1),modelts,'UniformOutput',false);
elseif catiter
    datats  = cellfun(@(C) cellcatrow(C,3),datats,'UniformOutput',false);
    modelts = cellfun(@(C) cellcatrow(C,3),modelts,'UniformOutput',false);
end
if catchan
  try
    if catrun && catiter
      datats  = cat(2,datats{:});
      modelts = cat(2,modelts{:});
    else
      datats  = cat(3,datats{:});
      modelts = cat(3,modelts{:});
      for pp=1:numruns
      for itr=1:numitr
          datats{itr,pp,1}  = cat(2,datats{itr,pp,:});
          modelts{itr,pp,1} = cat(2,modelts{itr,pp,:});
      end
      end
      datats  = datats(:,:,1);
      modelts = modelts(:,:,1);
    end
  catch ex
    warning(ex.message);
    warning('Skip cocatenating outputs across channels');
  end
end

%%% sub functions
%-- Concatenate cell array across columns
function M = cellcatcolumn(I,dim)
    nrow = size(I,1);
    M = cell(nrow,1);
    for ii=1:nrow
        M{ii,1} = cat(dim,I{ii,:});
    end

%-- Concatenate cell array across rows
function M = cellcatrow(I,dim)
    nclm = size(I,2);
    M = cell(1,nclm);
    for ii=1:nclm
        M{1,ii} = cat(dim,I{:,ii});
    end

%-- Check if option is 'yes'
function f = intprtopt(opt)
    if islogical(opt)
        f = opt;
    else
        f = strcmpi(opt,'yes');
    end

