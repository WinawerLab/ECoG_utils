function figlist = ecog_plotGridSC(dat,whichGrid, opt)

% p = ecog_plotGridSC(data, whichGrid, options)
% Plot in the entire grid.
% 
% specs should include 'timelabel','channels'

% Dependency: plotGridSCCommon, SetDefault, cellstrfind

% 20190726 Yuasa
% 20201117 Yuasa: enable to apply mask (opt.plot.AlphaData)
% 20210209 Yuasa: change option name (YLim, YTick, YTickLabel -> CLim, CTick, CTickLabel)
% 20220222 Yuasa: use plotGridSCCommon

%% Set options
narginchk(2,inf);
SetDefault('opt.timelabel',arrayfun(@(x) sprintf('Grid %02d',x),size(dat,3),'UniformOutput',false));
SetDefault('opt.plot.colorMap','parula');
SetDefault('opt.plot.colorBar',true);
SetDefault('opt.plot.nSubPlots',[]);
SetDefault('opt.plot.fontSize',12);
if isfield(opt.plot,'YLim') && ~isfield(opt.plot,'CLim')
    opt.plot.CLim = opt.plot.YLim;              % for compatible
end
if isfield(opt.plot,'YTick') && ~isfield(opt.plot,'CTick')
    opt.plot.CTick = opt.plot.YTick;            % for compatible
end
if isfield(opt.plot,'YTickLabel') && ~isfield(opt.plot,'CTickLabel')
    opt.plot.CTickLabel = opt.plot.YTickLabel;  % for compatible
end
SetDefault('opt.plot.CLim','maxmin');
SetDefault('opt.plot.CTick','auto',true);
SetDefault('opt.plot.CTickLabel','auto',true,'cell');
SetDefault('opt.plot.FigName','');
SetDefault('opt.plot.RotGrid',true);
SetDefault('opt.plot.ShowLabel',true);

SetDefault('opt.plot.AlphaData',[]);
SetDefault('opt.plot.Properties',{},'cell');

% whichHDgrid     = upper(whichHDgrid);
%%
%-- check mask
if ~isempty(opt.plot.AlphaData)
  if ndims(opt.plot.AlphaData)==1
    assert(size(dat,1)==size(opt.plot.AlphaData,1),'Mask data must consist of the same channels as the data.');
    maskdat = repmat(opt.plot.AlphaData,[1,size(dat,[2:ndims(dat)])]);
  else
    assert(isequal(size(dat),size(opt.plot.AlphaData)),'Mask data must be the same size as the data.');
    maskdat = opt.plot.AlphaData;
  end
else
  maskdat = ones(size(dat));
end

%%% run plotGridSCCommon
plotGridSCCommon;

%-- Plot figures
for itim = 1:length(opt.timelabel)
  for ee = 1:length(inx)
    %-- Decide how many subplots are needed
    nRow = opt.plot.nSubPlots(1);
    nCol = opt.plot.nSubPlots(2);
    
    plotElectrodes = gridList(inx{ee});
    
    %%% construct image data & alpha data for NaN
    scidx = reshape(inx{ee},nCol,nRow)';
    sensidx = grididx(scidx);
    imdat = nan(nRow,nCol);
    imdat(~isnan(sensidx))   = dat(sensidx(~isnan(sensidx)),itim);
    imalp = ones(nRow,nCol);
    imalp(~isnan(sensidx))   = maskdat(sensidx(~isnan(sensidx)),itim);
    imalp(isnan(imdat))   = 0;
    
    %%% imsc plot
    hF = figure('Name', opt.plot.FigName); 
    set(gcf, 'Position', [300 400 round(nCol*50) round(nRow*50)]);
    him = imagesc(imdat); hsc = him.Parent; box off;
    colormap(hsc,opt.plot.colorMap);
    him.AlphaData = imalp;   clim(hsc,[yrange]);
    if opt.plot.colorBar,  hcb = colorbar('FontSize', opt.plot.fontSize);    end 
    if setytick,        set(hcb,'Ticks',opt.plot.CTick);    end
    if setyticklabel,   set(hcb,'TickLabels',opt.plot.CTickLabel);    end
    set(hsc,'Position',get(hsc,'Position').*[0 1 0 1]+[0.078+0.032*~opt.plot.colorBar 0 0.78 0],...
        'YAxisLocation','right','FontSize', opt.plot.fontSize,opt.plot.Properties{:});
      hbg = axes('xlim',0.5 + [0 size(imdat,2)],'ylim',0.5 + [0 size(imdat,1)],...
          'color','none','XAxisLocation','top','YAxisLocation','left',...
          'FontSize', opt.plot.fontSize,opt.plot.Properties{:});
      set(hbg,'Position',get(hsc,'Position'));  axis(hbg,'ij');
    
    if ~opt.plot.ShowLabel
        set(hbg,'XTick',[],'YTick',[]);
        set(hsc,'XTick',[],'YTick',[]);
    elseif opt.plot.RotGrid
        set(hbg,'YTickLabel',plotElectrodes(1:nCol:end),'XTick',[]);
        set(hsc,'YTickLabel',plotElectrodes(nCol:nCol:end),'XTick',[]);
    else
        set(hbg,'XTickLabel',plotElectrodes(1:nCol),'YTick',[]);
        set(hsc,'XTickLabel',plotElectrodes(end+(1-nCol:0)),'YTick',[]);
    end
    title(opt.timelabel(itim))
    figlist{ee,itim} = hF;
  end
end
