function figlist = ecog_plotGridPRFts(data,stimulus,result,whichHDgrid, opt)

% p = ecog_plotGridPRFts(data,stimulus,result,whichHDgrid, option, data1)
% p = ecog_plotGridPRFts(data,stimulus,result,whichElectrodes, option, data1)
% ECOG_PLOTGRIDPRFTS plot time series in the entire HD grid.
% 
% <whichHDgrid> can be 'GA' or 'GB'.
% Otherwise it wll be interpreted as electrodes names (some options will be ignored).
% 
% <data> is a cell-array of time-series data.
% <stimulus> is a cell-array of stimulus data.
% <result> is a structure of pRF information (output of analyzePRF or analyzePRFdog), requiring: 
%   result.ecc:    Nx1 array with pRF eccentricity
%   result.ang:    Nx1 array with pRF angle [degree: 0deg is on positive x-axis and increase anti-clockwise]
%   result.rfsize: Nx1 array with pRF size
% 
% <option> is a structure of plot information.
%   option.channels     = a table of channel information.
%  	option.viselec      = a structure of visual area information.
%   option.plot         = a structure of plot information with following fields.
%       pix2deg         = degree/pixel (default = 1)
%       XLim            =
%       YLim            =
%       colors          =
% 
%       boundaries      = 'yes', 'no'(default) or numels  % show boundaries of events
%       addR2ToTitle    = 'yes' or 'no'(default)  % show R2 or xval in result
%       addEccToTitle   = 'yes' or 'no'(default)
%       nSubPlots       = [m n], plot m-by-n grid electrodes in a figure.
%       RotGrid         = true or false (default), if true, electrodes are
%                         aligned along row instead of column.
%       labelCol        = color or cell-array of colors, color of the titles
%                         in each axis. defualt is black.
%       FigName         = string, name of the figure.
%       fontSize        = scalar, font size in the figure.
%       legend          = string or cell-array of strings, legend labels.
%       options         = any other options for each axis.
% 
% Following options refer to result.options, but you can specify in <option>
%   option.hrf
%   option.maxpolydeg
%   option.gaussianmode
%   option.targetBAND
% 
% 

% option should include 'channels'

% Dependency: drawellipse, SetDefault, cellstrfind, arrangeinrect

% 20200310 Yuasa

%% Set options
narginchk(4,inf);
if isfield(opt,'channels')&&~isempty(opt.channels)
    channels = opt.channels;
elseif isfield(result,'channels')&&~isempty(result.channels)
    channels = result.channels;
else
    error('data or option must have ''channels'' field');
end
SetDefault('opt.plot.addR2ToTitle', 'no', 0);
SetDefault('opt.plot.addEccToTitle', 'no', 0);
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
SetDefault('result.options.hrf',1,0)
SetDefault('opt.hrf',result.options.hrf,0);
SetDefault('result.options.maxpolydeg',1,0)
SetDefault('opt.maxpolydeg',result.options.maxpolydeg,0);
SetDefault('result.options.gaussianmode','og',0);
SetDefault('opt.gaussianmode',result.options.gaussianmode,0);
SetDefault('result.options.targetBAND','bb',0)
SetDefault('result.targetBAND',result.options.targetBAND,0)
SetDefault('opt.targetBAND',result.targetBAND,0);
SetDefault('result.options.smoothingMode','none',0)
SetDefault('result.smoothingMode',result.options.smoothingMode,0)
SetDefault('opt.smoothingMode',result.smoothingMode,0);
SetDefault('result.options.smoothingN',1,0)
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
res          = size(stimulus{1},[1,2]);
resmx        = max(res);
hrf          = opt.hrf;
degs         = opt.maxpolydeg;
gaussianmode = lower(opt.gaussianmode);
tarBAND      = opt.targetBAND;
switch opt.smoothingMode
    case {'decimate'},  Dsample = opt.smoothingN;
    otherwise,          Dsample = 1;
end

%-- plot color properties
if isempty(opt.plot.colors)
    plotCol = get(groot,'DefaultAxesColorOrder');
elseif isnumeric(opt.plot.colors)
    if (size(opt.plot.colors,2)~=3 || any(opt.plot.colors(:)>1))
        %-- specify as integer
        plotCol = zeros(numel(opt.plot.colors),3);
        ColList = get(groot,'DefaultAxesColorOrder');
        ColList = ColList([end,1:end-1],:);
        for icol=1:numel(opt.plot.colors)
            if opt.plot.colors(icol) == 0
                plotCol(icol,:) = [1 1 1].*0.3;
            else
                plotCol(icol,:) = ColList(mod(opt.plot.colors(icol),size(ColList,1))+1,:);
            end
        end
    else
        %-- specify as rgb
        plotCol = opt.plot.colors;
    end
elseif ischar(opt.plot.colors) ||  isstring(opt.plot.colors)
    %-- specify as 'r' 'g' 'b' or "#0072BD" "#D95319"
    plotCol = reshape(opt.plot.colors,[],1);
end 

%-- check HD grid
if ischar(whichHDgrid) && ismember(whichHDgrid,{'GA','GB'})
    hasHDgrid   = true;
    gridList    = [];
else
    hasHDgrid   = false;
    gridList    = channels.name(cellstrfind(channels.name, whichHDgrid,1));
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
    elseif ee > numel(gridList)
      gridList{ee} = 'none';
    end
    
    %-- get Visual Area Label
    igrid = find(strcmp(channels.name,gridList{ee}),1);
    if strcmpi(opt.plot.addR2ToTitle,'yes') && ~isempty(igrid)
        areaList{ee} = sprintf('%s (%.1f%%)\n',gridList{ee},result.(R2field)(igrid));
    else
        areaList{ee} = sprintf('%s\n',gridList{ee});
    end
    if ~isempty(igrid)
        VA = [];
        if isTableCol(channels,'wangarea')
           VA = channels.wangarea{igrid};
        elseif isfield(opt,'viselec') && isfield(opt.viselec,'wang15_mplbl')
           vigrid = find(strcmp(opt.viselec.wang15_mplbl.elec_labels,gridList{ee}),1);
           if ~isempty(vigrid),  VA = opt.viselec.wang15_mplbl.area_labels{vigrid};  end
        elseif isfield(opt,'viselec') && isfield(opt.viselec,'wang2015_atlas')
           vigrid = find(strcmp(opt.viselec.wang2015_atlas.elec_labels,gridList{ee}),1);
           if ~isempty(vigrid),  VA = opt.viselec.wang2015_atlas.area_labels{vigrid};  end
        end
        if ~isempty(VA)&&~strcmp(VA,'none'), areaList{ee} = sprintf('%s w:%s',areaList{ee},VA); end
        
        VA = [];
        if isTableCol(channels,'bensonarea')
           VA = channels.bensonarea{igrid};
        elseif isfield(opt,'viselec') && isfield(opt.viselec,'benson14_varea')
           vigrid = find(strcmp(opt.viselec.benson14_varea.elec_labels,gridList{ee}),1);
           if ~isempty(vigrid),  VA = opt.viselec.benson14_varea.area_labels{vigrid};  end
        end
        if ~isempty(VA)&&~strcmp(VA,'none')
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

%%% Model computation
%-- Pre-compute cache for faster execution
[~,xx,yy] = makegaussian2d(resmx,2,2,2,2);

%-- Prepare the stimuli for use in the model
stimulusPP = {};
for pp=1:length(stimulus)
  stimulusPP{pp} = squish(stimulus{pp},2)';  % this flattens the image so that the dimensionality is now frames x pixels
  stimulusPP{pp} = [stimulusPP{pp} pp*ones(size(stimulusPP{pp},1),1)];  % this adds a dummy column to indicate run breaks
end

%-- Set model
switch gaussianmode
    case {'dog','gs'}
        modelfun = @(pp,dd) conv2run(modeldogcss(pp(1:5),pp(6:end),dd,res,xx,yy,0,0),hrf,dd(:,prod(res)+1));
    case {'og'}
        modelfun = @(pp,dd) conv2run(posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res),makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*abs(pp(3))^2))); 0]) .^ posrect(pp(5)),hrf,dd(:,prod(res)+1));
end

%-- Construct projection matrices
polymatrix = {};
for pp=1:length(degs)
  polymatrix{pp} = projectionmatrix(constructpolynomialmatrix(size(data{pp},2),0:degs(pp)));
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
      plotElectrodes = gridList(inx{ee}(ichidx));
      nameElectrodes = areaList(inx{ee}(ichidx));
      el = cellstrfind(channels.name,plotElectrodes,0);
      if ~isempty(el)
          
          
          datats = {};
          modelts = {};
          for pp=1:length(data)
              datats{pp}  = negfit.*polymatrix{pp}*data{pp}(el,:)';
              modelts{pp} = negfit.*polymatrix{pp}*modelfun(result.params(1,:,el),stimulusPP{pp});
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
          hplts(2) = plot(cat(1,modelts{:}),'r-','LineWidth', 2, opt.plot.data_options{:});
          pl = straightline(boundaries,'v','g-');  set(pl,'LineWidth',1);
          
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
