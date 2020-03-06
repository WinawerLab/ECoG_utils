function trials_out = ecog_plotGrid(whichHDgrid, opt, varargin)

% trials_out = ECOG_PLOTGRID(whichHDgrid, option, plot_function, x1, x2, ...)
% trials_out = ECOG_PLOTGRID(whichHDgrid, option, plot_function1, x1, x2,..., plot_function2, y1, y2,...)
% 
% ECOG_PLOTGRID is general function to plot data using any plot_function in the entire HD grid
% x1, x2,... are Nx1 cell-array with arguments of plot_function 
% 
% option is a structure of plot information.
%   option.channels     = a table of channel information.
%  	option.viselec      = (option) a structure of visual area information.
%   option.plot         = (option) a structure of plot information with following fields.
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
% Example:
%     opt = [];
%     opt.channels     = channels;
%     opt.plot.RotGrid = true;
%     opt.plot.options = {'XTick',[1 10 100],'XTickLabel',int2str([1 10 100]')};
%     ecog_plotGrid('GB',opt,@loglog,f,spectra,'Color','#0072BD');

% Dependency: SetDefault, cellstrfind

% 20191112 Yuasa
% 20191210 Yuasa: minor update

%% Set options
narginchk(3,inf);
assert(isfield(opt,'channels')&&~isempty(opt.channels),'options.channels is required');

SetDefault('opt.plot.addEccToTitle', 'no', 0);
SetDefault('opt.plot.fontSize', 12, 0);
SetDefault('opt.plot.nSubPlots', [], 1);

SetDefault('opt.plot.FigName', '', 0);
SetDefault('opt.plot.RotGrid', false, 0);
SetDefault('opt.plot.labelCol', [], 1);
SetDefault('opt.plot.legend', {}, 1);
SetDefault('opt.plot.options', {}, 1);

whichHDgrid     = upper(whichHDgrid);
%%
%-- process arguments
nchan = size(opt.channels,1);
multidat = find(cellfun(@(x) isa(x,'function_handle'), varargin));
ndat = numel(multidat);
func = varargin(multidat);
args = cell(1,ndat);
for idat = 1:ndat
    if idat == ndat
        argidx = (multidat(idat)+1):length(varargin);
    else
        argidx = (multidat(idat)+1):(multidat(idat+1)-1);
    end
    for ivar = argidx
        if ~iscell(varargin{ivar})
            varargin{ivar} = varargin(ivar);
        end
        if numel(varargin{ivar})==1
            varargin{ivar} = repmat(varargin{ivar},nchan,1);
        end
        varargin{ivar} = reshape(varargin{ivar}(:),[],1);
    end
    args{idat} = cat(2,varargin{argidx});
end

%-- set grid parameters
switch whichHDgrid
    case 'GA',  nCol = 8; nRow = 8;
    case 'GB',  nCol = 8; nRow = 16;
end
if opt.plot.RotGrid, tmp = nCol; nCol = nRow; nRow = tmp; end
FullGRID = nCol * nRow;

if isempty(opt.plot.nSubPlots)
    nFig = ceil(nRow / nCol ./ 1.2);
    opt.plot.nSubPlots = [ceil(nRow./nFig), nCol];
end
plRow  = opt.plot.nSubPlots(1);
nFig = ceil(nRow./plRow);

%-- correct elecnames (set '%03d')
[eleccat, elecnum] = strtok(opt.channels.name,int2str(0:9));
for el = 1:length(opt.channels.name)
    opt.channels.name{el} = sprintf('%s%03d',eleccat{el},str2double(elecnum{el}));
end

%-- Create a list of electrode names for the overall grid layout
gridList = [];
areaList = [];
isTableCol = @(t, thisCol) ismember(thisCol, t.Properties.VariableNames);
for ee = 1:FullGRID
    gridList{ee} = sprintf('%s%03d',whichHDgrid,ee);
    
    %-- get Visual Area Label
    areaList{ee} = sprintf('%s\n',gridList{ee});
    igrid = find(strcmp(opt.channels.name,gridList{ee}),1);
    if ~isempty(igrid)
        VA = [];
        if isTableCol(opt.channels,'wangarea')
           VA = opt.channels.wangarea{igrid};
        elseif isfield(opt,'viselec') && isfield(opt.viselec,'wang15_mplbl')
           vigrid = find(strcmp(opt.viselec.wang15_mplbl.elec_labels,gridList{ee}),1);
           if ~isempty(vigrid),  VA = opt.viselec.wang15_mplbl.area_labels{vigrid};  end
        elseif isfield(opt,'viselec') && isfield(opt.viselec,'wang2015_atlas')
           vigrid = find(strcmp(opt.viselec.wang2015_atlas.elec_labels,gridList{ee}),1);
           if ~isempty(vigrid),  VA = opt.viselec.wang2015_atlas.area_labels{vigrid};  end
        end
        if ~isempty(VA)&&~strcmp(VA,'none'), areaList{ee} = sprintf('%s w:%s',areaList{ee},VA); end
        
        VA = [];
        if isTableCol(opt.channels,'bensonarea')
           VA = opt.channels.bensonarea{igrid};
        elseif isfield(opt,'viselec') && isfield(opt.viselec,'benson14_varea')
           vigrid = find(strcmp(opt.viselec.benson14_varea.elec_labels,gridList{ee}),1);
           if ~isempty(vigrid),  VA = opt.viselec.benson14_varea.area_labels{vigrid};  end
        end
        if ~isempty(VA)&&~strcmp(VA,'none')
            areaList{ee} = sprintf('%s b:%s',areaList{ee},VA);
            switch opt.plot.addEccToTitle
                case 'yes'
                    if isTableCol(opt.channels,'bensonarea')
                        VC = opt.channels.bensoneccen(igrid);
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
if ~opt.plot.RotGrid
    v = reshape(1:FullGRID,[nRow,nCol]);
else
    v = flipud(reshape(1:FullGRID,[nCol,nRow])'); 
end
inx = [];
for ee = 1:nFig
    inx{ee} = reshape(v(plRow*(ee-1)+1:min(nRow,plRow*(ee)),:)',1,[]);
end

%-- legend
lgnds = opt.plot.legend;
if ~iscell(lgnds), lgnds = {lgnds}; end
flglgnd = ~isempty(lgnds);

%-- main
for ee = 1:length(inx)
    %-- Decide how many subplots are needed
    if ~isempty(opt.plot.nSubPlots)
        nRow = opt.plot.nSubPlots(1);
        nCol = opt.plot.nSubPlots(2);
    else
        nPlot = length(el_index);
        nRow = ceil(sqrt(nPlot));
        nCol = ceil(sqrt(nPlot));
        if nPlot <= (nRow*nCol)-nCol
            nRow = nRow-1;
        end
    end
       
    %-- subplot
    hF = figure('Name', opt.plot.FigName); 
    set(hF, 'Position', [150 100 2000 1250]);
    for ichidx=1:length(inx{ee})
      plotElectrodes = gridList(inx{ee}(ichidx));
      nameElectrodes = areaList(inx{ee}(ichidx));
      el = cellstrfind(opt.channels.name,plotElectrodes,0);
      if ~isempty(el)
                    
          %-- plot       hA = subplot(nRow,nCol,ichidx);
          set(0,'CurrentFigure',hF);
          ax_sp = [0.03 0.03]; ax_siz = (1-ax_sp*2) ./[nCol nRow];
          ax_pos = [ax_sp(1) + ax_siz(1)*(mod(ichidx-1,nCol)),...
                    ax_sp(2) + ax_siz(2)*(nRow-ceil(ichidx./nCol)),...
                    ax_siz(1)*0.70, ax_siz(2)*0.70];
          hA = subplot('Position',ax_pos);
              outerpos = hA.OuterPosition;
              ax_ti = hA.TightInset;
              ax_left = outerpos(1) + ax_ti(1);
              ax_bottom = outerpos(2) + ax_ti(2)*2;
              ax_width = outerpos(3) - ax_ti(1) - ax_ti(3);
              ax_height = outerpos(4) - ax_ti(2)*2 - ax_ti(4)*4;
              hA.Position = [ax_left ax_bottom ax_width ax_height];
              
          hplts = [];
          for idat = 1:ndat
              feval(func{idat},args{idat}{el,:});
              hold on;
              hplts(idat) = hA.Children(1);
          end 
          hold off;
          set(hA,'FontSize', opt.plot.fontSize,opt.plot.options{:});
          title(nameElectrodes);
          if flglgnd
              lgndsidx = ~cellfun(@isempty,lgnds);
              hlg = legend(hplts(lgndsidx),lgnds(lgndsidx),...
                  'Location','southwest','FontSize',opt.plot.fontSize);
              hlg.Position = hlg.Position .* [0 1 1 1] + [hA.OuterPosition(1)-hlg.Position(3)-0.01 0 0 0];
              flglgnd = false;
          end
          drawnow limitrate nocallbacks;
      end
    end
    trials_out{ee} = hF;
end
