function figlist = ecog_plotGridPRFts(data,stimulus,result,whichHDgrid, opt)

% p = ecog_plotGridPRFts(data,stimulus,result,whichHDgrid, option)
% p = ecog_plotGridPRFts(data,stimulus,result,whichElectrodes, option)
% p = ecog_plotGridPRFts(data,stimulus,params,whichHDgrid, option)
% p = ecog_plotGridPRFts(data,stimulus,params,whichElectrodes, option)
% ECOG_PLOTGRIDPRFTS plot time series in the entire HD grid.
% 
% <whichHDgrid> can be 'GA' or 'GB'.
% Otherwise it wll be interpreted as electrodes names (some options will be ignored).
% 
% <data> is a cell-array of time-series data.
% <stimulus> is a cell-array of stimulus data.
% <result> is a structure of pRF information (output of analyzePRF or analyzePRFdog), requiring: 
%   result.params:    1xNxP array with pRF parameters
%   result.options:   option structure to estimate pRFs
% <params> is result.params, which you can directly specify instead of result.
% 
% <option> is a structure of plot information.
%   option.channels      = a table of channel information.
%                          If result has channel field, it's not necessary.
%  	option.viselec       = a structure of visual area information.
%                          If channels has visual area table or you don't want to show area name,
%                          it's not necessary.
%   option.plot          = a structure of plot information with following fields.
%       pix2deg          = degree/pixel (default = 1) (required if option.plotbenson = 'yes')
%       XLim             =
%       YLim             =
%       colors           =
% 
%       boundaries       = 'yes', 'no'(default) or numels  % show boundaries of events
%       addChsToTitle    = 'yes'(default) or 'no'
%       addWangToTitle   = 'yes'(default) or 'no'
%       addBensonToTitle = 'yes'(default) or 'no'
%       addR2ToTitle     = 'yes' or 'no'(default)  % show R2 or xval of result
%       addEccToTitle    = 'yes' or 'no'(default)  % show eccentricity of result
%       addSbjToTile     = 'yes' or 'no'(default)  % show subject name
%       nSubPlots        = [m n], plot m-by-n grid electrodes in a figure.
%       RotGrid          = true or false (default), if true, electrodes are
%                          aligned along row instead of column.
%       labelCol         = color or cell-array of colors, color of the titles
%                          in each axis. defualt is black.
%       FigName          = string, name of the figure.
%       fontSize         = scalar, font size in the figure.
%       legend           = string or cell-array of strings, legend labels.
%       options          = any other options for each axis.
%   option.plotbenson    = 'yes' or 'no'(default)  % plot results estimated from benson atlas
%   option.subjet        = string, subject name is helpful to interpret benson prf.
% 
% Following options refer to result.options, but you can specify in <option>
%   option.hrf
%   option.maxpolydeg
%   option.gaussianmode
%   option.targetBAND
% 
% 

% Dependency: <analyzePRF>, SetDefault, cellstrfind, arrangeinrect, istablefield

% 20200310 Yuasa
% 20200331 Yuasa: add 'plotbenson' option
% 20200925 Yuasa: solve problem using strings in channels

%% Set options
narginchk(4,inf);
if isfield(opt,'channels')&&~isempty(opt.channels)
    channels = opt.channels;
elseif isfield(result,'channels')&&~isempty(result.channels)
    channels = result.channels;
else
    error('data or option must have ''channels'' field');
end
if ~isstruct(result),   result = struct('params',result);   end

SetDefault('opt.plotbenson', 'no', 0);

SetDefault('opt.plot.addChsToTitle', 'yes', 0);
SetDefault('opt.plot.addWangToTitle', 'yes', 0);
SetDefault('opt.plot.addBensonToTitle', 'yes', 0);
SetDefault('opt.plot.addR2ToTitle', 'no', 0);
SetDefault('opt.plot.addEccToTitle', 'no', 0);
SetDefault('opt.plot.addSbjToTitle', 'no', 0);
SetDefault('opt.plot.fontSize', 12, 0);
SetDefault('opt.plot.pix2deg', 1, 1);
SetDefault('opt.plot.nSubPlots', [], 1);
SetDefault('opt.plot.XLim', [], 1);
SetDefault('opt.plot.YLim', [], 1);
SetDefault('opt.plot.boundaries', 'no', 0);

SetDefault('opt.plot.FigName', '', 0);
SetDefault('opt.plot.RotGrid', false, 0);
SetDefault('opt.plot.colors', [], 1);
SetDefault('opt.plot.labelCol', [], 1);
SetDefault('opt.plot.representative', [], 1);
SetDefault('opt.plot.legend', {}, 1,'cell');
SetDefault('opt.plot.data_options', {}, 1);
SetDefault('opt.plot.model_options', {}, 1);

SetDefault('whichHDgrid','*',0);

%-- Set model parameters
SetDefault('result.options.hrf',1,0);
SetDefault('opt.hrf',result.options.hrf,0);
SetDefault('result.options.maxpolydeg',1,0);
SetDefault('opt.maxpolydeg',result.options.maxpolydeg,0);
SetDefault('result.options.gaussianmode','og',0);
SetDefault('opt.gaussianmode',result.options.gaussianmode,0);
SetDefault('result.options.targetBAND','bb',0);
SetDefault('result.targetBAND',result.options.targetBAND,0)
SetDefault('opt.targetBAND',result.targetBAND,0);
SetDefault('result.options.smoothingMode','none',0);
SetDefault('result.smoothingMode',result.options.smoothingMode,0)
SetDefault('opt.smoothingMode',result.smoothingMode,0);
SetDefault('result.options.smoothingN',1,0);
SetDefault('result.smoothingN',result.options.smoothingN,0)
SetDefault('opt.smoothingN',result.smoothingN,0);

%-- Check boundaries
if ischar(opt.plot.boundaries) && strcmp(opt.plot.boundaries,'yes')
    assert(isfield(opt,'events'),'option.plot.boundaries = ''yes'' requires option.events');
end

%-- convert into cell
if ~iscell(data), data = {data}; end
if ~iscell(stimulus), stimulus = {stimulus}; end

%%
%-- Define variables
ndats        = 1;
hrf          = opt.hrf;
degs         = opt.maxpolydeg;
gaussianmode = lower(opt.gaussianmode);
tarBAND      = opt.targetBAND;
switch opt.smoothingMode
    case {'decimate'},  Dsample = opt.smoothingN;
    otherwise,          Dsample = 1;
end
plotbenson   = false;
if strcmpi(opt.plotbenson,'yes')
    if all(ismember({'bensoneccen','bensonangle','bensonsigma'},channels.Properties.VariableNames))
        plotbenson   = true;
        ndats = ndats + 1;
    else
        warning('Benson atlas information is not found');
    end
end
if ~istablefield(channels,'subject_name')
    if isfield(opt,'subject'),        channels.subject_name(:) = opt.subject;
    else
        if isfield(result,'subject'), channels.subject_name(:) = result.subject;
        else,                         channels.subject_name(:) = '';
        end
    end
end
subjects = cellstr(channels.subject_name);
numvxs = height(channels);
numruns = size(data,2);

%-- Get data length
ntime = zeros(numvxs,numruns);
for pp=1:numruns
    ntime(:,pp) = sum(~isnan(data{1,pp}),2);
end

%-- Set stimulus (duplicate)
if ~iscell(stimulus{1}),   stimulus = {stimulus};  end
if length(stimulus)==1,    stimulus = repmat(stimulus,numvxs,1);    end

res = sizefull(stimulus{1}{1},2);
resmx = max(res);

%-- plot color properties
ColList = ['rbgymc']';
if isempty(opt.plot.colors)
    opt.plot.colors = reshape(1:size(ColList,1),[],1);
end
if isnumeric(opt.plot.colors)
    if (size(opt.plot.colors,2)~=3 || any(opt.plot.colors(:)>1))
        %-- specify as integer
        plotCol = ColList(mod(opt.plot.colors(:)-1,size(ColList,1))+1,:);
        plotCol(opt.plot.colors==0,:) = 'k';   % special color for '0'
    else
        %-- specify as rgb
        plotCol = opt.plot.colors;
    end
elseif ischar(opt.plot.colors) ||  isstring(opt.plot.colors)
    %-- specify as 'r' 'g' 'b' or "#0072BD" "#D95319"
    plotCol = reshape(opt.plot.colors,[],1);
end 
plotCol = repmat(plotCol,ceil(ndats./size(plotCol,1)),1);

%-- check HD grid
if ischar(whichHDgrid) && ismember(whichHDgrid,{'GA','GB'})
    hasHDgrid   = true;
    gridList    = [];
    chanIdx     = [];
else
    hasHDgrid   = false;
    if islogical(whichHDgrid)
      chanIdx     = find(whichHDgrid);
    elseif isnumeric(whichHDgrid) && isequaln(whichHDgrid,round(whichHDgrid))
      chanIdx     = whichHDgrid;
    else    % list of channel names
      chanIdx     = cellstrfind(channels.name, whichHDgrid,0);
    end
    gridList    = channels.name(chanIdx);
    whichHDgrid = 'Channels';
end

%-- correct elecnames (set '%03d')
if hasHDgrid
[eleccat, elecnum] = strtok(channels.name,int2str(0:9));
for el = 1:length(channels.name)
    channels.name{el} = sprintf('%s%03d',eleccat{el},str2double(elecnum{el}));
end
end

%-- set grid parameters
switch whichHDgrid
    case 'GA',  nCol = 8; nRow = 8;
    case 'GB',  nCol = 8; nRow = 16;
    otherwise,  [nCol,nRow] = arrangeinrect(numel(gridList),1.0,[0.6,2.0]);
end
if hasHDgrid && opt.plot.RotGrid, tmp = nCol; nCol = nRow; nRow = tmp; end
FullGRID = nCol * nRow;

if isempty(opt.plot.nSubPlots)
    nFig = ceil(nRow / nCol ./ 1.7);
    opt.plot.nSubPlots = [ceil(nRow./nFig), nCol];
end
plRow  = opt.plot.nSubPlots(1);
nFig = ceil(nRow./plRow);

%-- Check cross-validation
if isfield(result,'xval') && ~all(isnan(result.xval))
    R2field = 'xval';
elseif isfield(result,'aggregatedtestperformance') && ~all(isnan(result.aggregatedtestperformance))
    R2field = 'aggregatedtestperformance';
else
    R2field = 'R2';
end

%-- Create a list of electrode names for the overall grid layout
areaList = [];
isTableCol = @(t, thisCol) ismember(thisCol, t.Properties.VariableNames);
for ee = 1:FullGRID
    if hasHDgrid
      gridList{ee} = sprintf('%s%03d',whichHDgrid,ee);
      igrid = find(strcmp(channels.name,gridList{ee}),1);
      if isempty(igrid),    chanIdx(ee)  = nan;
      else,                 chanIdx(ee)  = igrid;
      end
    elseif ee > numel(gridList)
      gridList{ee} = 'none';
      chanIdx(ee)  = nan;
    end
    
    %-- get Visual Area Label
    igrid = chanIdx(ee);
    if strcmpi(opt.plot.addR2ToTitle,'yes') && ~isnan(igrid)
        if strcmpi(opt.plot.addChsToTitle,'yes')
            areaList{ee} = sprintf('%s (%.1f%%)',gridList{ee},result.(R2field)(igrid));
        else
            areaList{ee} = sprintf('(%.1f%%)',gridList{ee},result.(R2field)(igrid));
        end
    elseif strcmpi(opt.plot.addChsToTitle,'yes')
        areaList{ee} = sprintf('%s',gridList{ee});
    else
        areaList{ee} = '';
    end
    if (strcmpi(opt.plot.addWangToTitle,'yes')||strcmpi(opt.plot.addBensonToTitle,'yes'))&& ...
            (strcmpi(opt.plot.addChsToTitle,'yes')||strcmpi(opt.plot.addR2ToTitle,'yes'))
        areaList{ee} = sprintf('%s\n',areaList{ee});
    end
    if ~isnan(igrid)
        if strcmpi(opt.plot.addSbjToTitle,'yes') && ~isempty(subjects{ee})
            areaList{ee} = sprintf('%s %s',subjects{igrid}, areaList{ee});
        end
        
        VA = [];
        if isTableCol(channels,'wangarea')
           if iscell(channels.wangarea),    VA = channels.wangarea{igrid};
           else,                            VA = channels.wangarea(igrid);
           end
        elseif isfield(opt,'viselec') && isfield(opt.viselec,'wang15_mplbl')
           vigrid = find(strcmp(opt.viselec.wang15_mplbl.elec_labels,gridList{ee}),1);
           if ~isempty(vigrid),  VA = opt.viselec.wang15_mplbl.area_labels{vigrid};  end
        elseif isfield(opt,'viselec') && isfield(opt.viselec,'wang2015_atlas')
           vigrid = find(strcmp(opt.viselec.wang2015_atlas.elec_labels,gridList{ee}),1);
           if ~isempty(vigrid),  VA = opt.viselec.wang2015_atlas.area_labels{vigrid};  end
        end
        if strcmpi(opt.plot.addWangToTitle,'yes')&&~isempty(VA)&&~strcmp(VA,'none')
            areaList{ee} = sprintf('%s w:%s',areaList{ee},VA);
        end
        
        VA = [];
        if isTableCol(channels,'bensonarea')
           if iscell(channels.bensonarea),  VA = channels.bensonarea{igrid};
           else,                            VA = channels.bensonarea(igrid);
           end
        elseif isfield(opt,'viselec') && isfield(opt.viselec,'benson14_varea')
           vigrid = find(strcmp(opt.viselec.benson14_varea.elec_labels,gridList{ee}),1);
           if ~isempty(vigrid),  VA = opt.viselec.benson14_varea.area_labels{vigrid};  end
        end
        if strcmpi(opt.plot.addBensonToTitle,'yes')&&~isempty(VA)&&~strcmp(VA,'none')
            areaList{ee} = sprintf('%s b:%s',areaList{ee},VA);
            switch opt.plot.addEccToTitle
                case 'yes'
                    if isTableCol(channels,'bensoneccen')
                        VC = channels.bensoneccen(igrid);
                    elseif isfield(opt,'viselec') && isfield(opt.viselec,'benson14_varea')
                        VC = opt.viselec.benson14_varea.node_eccen(vigrid);
                    end
                    areaList{ee} = sprintf('%s [ecc = %0.1f]',areaList{ee},VC);
            end
        end        
    end
    labelCol = 'black';
    if ~isempty(opt.plot.labelCol)
        if size(opt.plot.labelCol,1)==1
            if iscell(opt.plot.labelCol), labelCol = opt.plot.labelCol{1};
            else,                         labelCol = opt.plot.labelCol;
            end
        else
            if iscell(opt.plot.labelCol), labelCol = opt.plot.labelCol{igrid};
            else,                         labelCol = opt.plot.labelCol(igrid,:);
            end
        end
    end
    areaList{ee} = sprintf('{\\color{%s} \\bf %s}',labelCol,areaList{ee});
end
 
%-- Plot response per electrode, in two separate figures (top and bottom)
if ~hasHDgrid
    v = reshape(1:FullGRID,[nCol,nRow])';
elseif ~opt.plot.RotGrid
    v = reshape(1:FullGRID,[nRow,nCol]);
else
    v = flipud(reshape(1:FullGRID,[nCol,nRow])'); 
end
inx = [];
for ee = 1:nFig
    inx{ee} = reshape(v(plRow*(ee-1)+1:min(nRow,plRow*(ee)),:)',1,[]);
end

%%% Benson reference
if plotbenson
    PRF2CSS  = @(pp,res,gain,expt) [(1+res)/2 - pp(1)*sin(pp(2)*pi/180),...
               pp(1)*sin(pp(2)*pi/180) + (1+res)/2,...
               pp(3)*expt, gain, expt, 1, 0];

    subj_L = {'beilen','som648','som692','som708','som748'}; % left hemisphere
    subj_R = {'chaam','som661','som674','som718','som726'}; % right hemisphere
    subj_B = {'som723'}; % bilateral
    assert(~isempty(setdiff(subjects,[subj_L,subj_R,subj_B])),...
            'data includes unknown subject');
        
    %-- detect channels in left hemisphere
    islefths = ismember(subjects,subj_L) | ...
                (ismember(subjects,subj_B) & startsWith(channels.name,'L'));
    bensonPRF = [channels.bensoneccen./opt.plot.pix2deg, ...
                 ((-1).^islefths).*channels.bensonangle+90, ...
                 channels.bensonsigma./opt.plot.pix2deg];

    bensonparams = zeros(1,7,size(channels,1));
    for el=1:size(channels,1)
%         bensonparams(:,:,el) = PRF2CSS(bensonPRF(el,:),resmx,result.params(1,4,el),result.params(1,5,el));
        bensonparams(:,:,el) = PRF2CSS(bensonPRF(el,:),resmx,1,1);    % gain = 1; expt = 1; (linear model with dummy gain)
    end
end

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
        modelfun = @(pp,dd) conv2run(posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res),makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*abs(pp(3))^2))); 0]) .^ posrect(pp(5)),hrf,dd(:,prod(res)+1));
end

%-- Construct projection matrices
polymatrix = repmat({{}},numvxs,1);
for vxs=1:numvxs
for pp=1:numruns
  polymatrix{vxs}{pp} = projectionmatrix(constructpolynomialmatrix(ntime(vxs,pp),0:degs(pp)));
end
end

%-- Make boundaries
if isnumeric(opt.plot.boundaries)
    boundaries = opt.plot.boundaries(:);
elseif ischar(opt.plot.boundaries) && strcmp(opt.plot.boundaries,'yes')
    stimtypelist = grp2idx(cellfun(@(x) strtok(x,'-'),opt.events.trial_name,'UniformOutput',false));
    boundaries = [0 reshape(find(diff(stimtypelist)~=0),1,[]) height(opt.events)] ./ Dsample;
else
    boundaries = [];
end

%-- reverse for alpha suppression
switch tarBAND
    case {'a','aC','aCb','aCR','aCRb','aCL','aCLb','FaC','FaCb','FaCR','FaCRb','FaCL','FaCLb'}
        negfit = -1;
        axmode = 'ij';
    otherwise
        negfit = 1;
        axmode = 'xy';
end


%%% main
for ee = 1:length(inx)
    %-- legend
    flglgnd = ~isempty(opt.plot.legend);
    %-- Decide how many subplots are needed
    nRow = opt.plot.nSubPlots(1);
    nCol = opt.plot.nSubPlots(2);
       
    %-- subplot
    hF = figure('Name', opt.plot.FigName); 
    set(hF, 'Position', [150 100 200*nCol 160*nRow]);%[150 100 2000 1250]
%             set(gcf,'Units','points');
    for ichidx=1:length(inx{ee})
%       plotElectrodes = gridList{inx{ee}(ichidx)};
      nameElectrodes = areaList{inx{ee}(ichidx)};
      el = chanIdx(inx{ee}(ichidx));
      if ~isnan(el)
          datats = {};
          modelts = {};
          for pp=1:numruns
              nanplace    = isnan(data{pp}(el,:));
              datats{pp}  = negfit.*polymatrix{el}{pp}*data{pp}(el,~nanplace)';
              modelts{pp} = negfit.*polymatrix{el}{pp}*modelfun(result.params(1,:,el),stimulusPP{el}{pp}(~nanplace,:));
          end
            if plotbenson
              modeltsB = {};
              for pp=1:numruns
                  modeltsB{pp} = negfit.*polymatrix{el}{pp}*modelfun(bensonparams(1,:,el),stimulusPP{el}{pp}(~nanplace,:));
                  bensonparams(1,4,el) = prctile(posrect(negfit.*datats{pp}),95)./prctile(posrect(negfit.*modeltsB{pp}),95);
                  modeltsB{pp} = negfit.*polymatrix{el}{pp}*modelfun(bensonparams(1,:,el),stimulusPP{el}{pp}(~nanplace,:));
              end
            end
          
          
          %-- Prepare axis       hA = subplot(nRow,nCol,ichidx);
          set(0,'CurrentFigure',hF);
          ax_sp = 1./([nCol nRow]).*[0.3 0.25];
          ax_siz = (1-ax_sp.*[1 1.2]) ./[nCol nRow];
          ax_pos = [ax_sp(1) + ax_siz(1)*(mod(ichidx-1,nCol)),...
                    ax_sp(2) + ax_siz(2)*(nRow-ceil(ichidx./nCol)),...
                    ax_siz(1)*0.80, ax_siz(2)*0.50];
          hA = subplot('Position',ax_pos);
          
          hplts = [];
          hplts(1) = plot(cat(1,datats{:}),'k:o', 'LineWidth', 2, 'MarkerFaceColor','auto', opt.plot.data_options{:});
          hold on;
          pl = straightline(0,'h','k:');  set(pl,'LineWidth',1);
          hplts(2) = plot(cat(1,modelts{:}),'-','Color',plotCol(1,:),'LineWidth', 2, opt.plot.data_options{:});
            if plotbenson
              hplts(3) = plot(cat(1,modeltsB{:}),'-','Color',plotCol(2,:),'LineWidth', 2, opt.plot.data_options{:});
            end
          pl = straightline(boundaries,'v','-');  set(pl,'Color',[1 1 1].*0.5,'LineWidth',1);
          
          %-- Set axis, title & legend
          if isempty(opt.plot.XLim),  xlims = [0 size(cat(1,modelts{:}),1)]+0.5;
          else,                       xlims = opt.plot.XLim;
          end
          if isempty(opt.plot.YLim),  ylims = ylim;
          else,                       ylims = opt.plot.YLim;
          end
          
          axis(hA,axmode);
          set(hA,'XLim',xlims,'YLim',ylims, 'FontSize', opt.plot.fontSize);
          title(nameElectrodes);
%           setsubplotaxes(hA);
          if flglgnd
              numlg = min(length(hplts),length(opt.plot.legend));
              hlg = legend(hplts(1:numlg),opt.plot.legend(1:numlg),'FontSize',opt.plot.fontSize);
              hlg.NumColumns = nCol;
              hlg.Position([1,2]) = [ax_sp(1), hA.OuterPosition([2,4])*[1;1.02]];
              flglgnd = false;
          end
          drawnow limitrate nocallbacks;
      end
    end
    figlist{ee} = hF;
end
