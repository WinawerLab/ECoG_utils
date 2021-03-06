function figlist = ecog_plotGridSC(dat,whichHDgrid, opt)

% p = ecog_plotGridSC(data, whichHDgrid, options)
% Plot the entire HD grid
% 
% specs should include 'timelable','channels','viselec'

% Dependency: SetDefault, cellstrfind

% 20190726 Yuasa
% 20201117 Yuasa - enable to apply mask (opt.plot.AlphaData)

%% Set options
narginchk(2,inf);
SetDefault('opt.plot.colorMap','parula');
SetDefault('opt.plot.colorBar',true);
SetDefault('opt.plot.nSubPlots',[]);
SetDefault('opt.plot.fontSize',12);
SetDefault('opt.plot.YLim','maxmin');
SetDefault('opt.plot.YTick','auto',true);
SetDefault('opt.plot.YTickLabel','auto',true,'cell');
SetDefault('opt.plot.FigName','');
SetDefault('opt.plot.RotGrid',true);

SetDefault('opt.plot.AlphaData',[]);
SetDefault('opt.plot.Properties',{},'cell');

% whichHDgrid     = upper(whichHDgrid);
%%
%-- check mask
setmask = ~isempty(opt.plot.AlphaData);
if setmask
    assert(isequal(size(dat),size(opt.plot.AlphaData)),'Mask data must be the same size as the data.');
    maskdat = opt.plot.AlphaData;
end
%-- set grid parameters
switch whichHDgrid
    case 'GA',  nCol = 8; nRow = 8;
    case 'GB',  nCol = 8; nRow = 16;
end
if opt.plot.RotGrid, tmp = nCol; nCol = nRow; nRow = tmp; end
FullGRID = nCol * nRow;

if isempty(opt.plot.nSubPlots)
    nFig = ceil(nRow / nCol ./ 2.2);
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
gridList = cell(1,FullGRID);
grididx  = nan(1,FullGRID);
for ee = 1:FullGRID
    gridList{ee} = sprintf('%s%03d',whichHDgrid,ee);
    idx = cellstrfind(opt.channels.name,gridList{ee});
    if ~isempty(idx), grididx(ee) = idx; end
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

%-- Set ylim
ichs = cellstrfind(opt.channels.name,gridList,1);
ymin = prctile(reshape(dat(ichs,:),1,[]),5);
ymax = prctile(reshape(dat(ichs,:),1,[]),95);
if ischar(opt.plot.YLim)
  switch opt.plot.YLim
    case 'maxmin',  yrange = [ymin ymax];
    case 'maxabs',  yrange = [-1 1].*max(abs([ymin ymax]));
    case 'minabs',  yrange = [-1 1].*min(abs([ymin ymax]));
    otherwise,      error('unknown parameter');
  end
else
    yrange = opt.plot.YLim;
end
setytick = opt.plot.colorBar && ...
    ~(ischar(opt.plot.YTick)&&strcmp(opt.plot.YTick,'auto'));
setyticklabel = opt.plot.colorBar && ...
    ~(ischar(opt.plot.YTickLabel{1})&&strcmp(opt.plot.YTickLabel{1},'auto'));

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
    if setmask
    imalp(~isnan(sensidx))   = maskdat(sensidx(~isnan(sensidx)),itim);
    end
    imalp(isnan(imdat))   = 0;
    
    %%% imsc plot
    hF = figure('Name', opt.plot.FigName); 
    set(gcf, 'Position', [300 400 round(nCol*50) round(nRow*50)]);
    him = imagesc(imdat); hsc = gca; box off;
    colormap(hsc,opt.plot.colorMap);
    him.AlphaData = imalp;   clim(hsc,[yrange]);
    if opt.plot.colorBar,  hcb = colorbar('FontSize', opt.plot.fontSize);    end 
    if setytick,        set(hcb,'Ticks',opt.plot.YTick);    end
    if setyticklabel,   set(hcb,'TickLabels',opt.plot.YTickLabel);    end
    set(hsc,'Position',get(hsc,'Position').*[0 1 0 1]+[0.078 0 0.78 0],...
        'YAxisLocation','right','FontSize', opt.plot.fontSize,opt.plot.Properties{:});
      hbg = axes('xlim',0.5 + [0 size(imdat,2)],'ylim',0.5 + [0 size(imdat,1)],...
          'color','none','XAxisLocation','top','YAxisLocation','left',...
          'FontSize', opt.plot.fontSize,opt.plot.Properties{:});
      set(hbg,'Position',get(hsc,'Position'));
    
    if opt.plot.RotGrid
        set(hbg,'YTickLabel',plotElectrodes(fliplr(1:nCol:end)),'XTick',[]);
        set(hsc,'YTickLabel',plotElectrodes(nCol:nCol:end),'XTick',[]);
    else
        set(hbg,'XTickLabel',plotElectrodes(1:nCol),'YTick',[]);
        set(hsc,'XTickLabel',plotElectrodes(end+(1-nCol:0)),'YTick',[]);
    end
    title(opt.timelabel(itim))
    figlist{ee,itim} = hF;
  end
end
