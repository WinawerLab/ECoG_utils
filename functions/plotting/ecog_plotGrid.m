function figlist = ecog_plotGrid(whichGrid, opt, varargin)

% p = ECOG_PLOTGRID(whichGrid, option, plot_function, x1, x2, ...)
% p = ECOG_PLOTGRID(whichGrid, option, plot_function1, x1, x2,..., plot_function2, y1, y2,...)
% 
% ECOG_PLOTGRID is general function to plot data using any plot_function in the entire grid
% x1, x2,... are Nx1 cell-array with arguments of plot_function 
% 
% option is a structure of plot information.
%   option.channels     = a table of channel information.
%  	option.viselec      = (option) a structure of visual area information.
%  	option.R2           = (option) an array of R^2 values.
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

% Dependency: plotGridCommon, SetDefault, cellstrfind

% 20191112 Yuasa
% 20191210 Yuasa: minor update
% 20210701 Yuasa: use plotGridCommon

%% Set options
narginchk(3,inf);
assert(isfield(opt,'channels')&&~isempty(opt.channels),'options.channels is required');

SetDefault('opt.plot.addChsToTitle', 'yes', 0);
SetDefault('opt.plot.addWangToTitle', 'yes', 0);
SetDefault('opt.plot.addBensonToTitle', 'yes', 0);
SetDefault('opt.plot.addR2ToTitle', 'no', 0);
SetDefault('opt.plot.addEccToTitle', 'no', 0);
SetDefault('opt.plot.addSbjToTitle', 'no', 0);
SetDefault('opt.plot.fontSize', 12, 0);
SetDefault('opt.plot.nSubPlots', [], 1);

SetDefault('opt.plot.FigName', '', 0);
SetDefault('opt.plot.RotGrid', false, 0);
SetDefault('opt.plot.labelCol', [], 1);
SetDefault('opt.plot.legend', {}, 1);
SetDefault('opt.plot.options', {}, 1);

% whichGrid     = upper(whichGrid);
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

channels = opt.channels;
subjects = cellstr(channels.subject_name);

%%% run plotGridCommon
plotGridCommon;

%-- legend
lgnds = opt.plot.legend;
if ~iscell(lgnds), lgnds = {lgnds}; end
flglgnd = ~isempty(lgnds);

%-- main
for ee = 1:length(inx)
    %-- Decide how many subplots are needed
    nRow = opt.plot.nSubPlots(1);
    nCol = opt.plot.nSubPlots(2);
       
    %-- subplot
    hF = figure('Name', opt.plot.FigName); 
    set(hF, 'Position', [150 100 2000 1250]);
    for ichidx=1:length(inx{ee})
      plotElectrodes = gridList(inx{ee}(ichidx));
      nameElectrodes = areaList(inx{ee}(ichidx));
      el = cellstrfind(channels.name,plotElectrodes,0);
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
    figlist{ee} = hF;
end
