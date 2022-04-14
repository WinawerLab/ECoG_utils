function figlist = ecog_plotGridXSCts(dat,whichGrid, opt)

% p = ecog_plotGridXSCts(data, whichGrid, options)
% Plot in the entire grid.
% 
% data is channels x channels x time matrix (considerd as multiple columns with each seed)
% specs requires following fields: 'timelabel','channels','seed'

% Dependency: plotGridSCCommon, SetDefault, cellstrfind

% 20210517 Yuasa
% 20220222 Yuasa: use plotGridSCCommon

%% Set options
narginchk(2,inf);
SetDefault('opt.timelabel',arrayfun(@(x) sprintf('Grid %02d',x),1:size(dat,3),'UniformOutput',false));
SetDefault('opt.seed',opt.channels.name{1});
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

SetDefault('opt.plot.Time',1);
SetDefault('opt.plot.ShowSeed',true);
SetDefault('opt.plot.Interaction',true);
SetDefault('opt.plot.AlphaData',[]);
SetDefault('opt.plot.Properties',{},'cell');

% whichHDgrid     = upper(whichHDgrid);
%%
%-- check mask
nchan = size(dat,1);
nxchan = size(dat,2);
assert(nxchan==nchan | nxchan==1,'1st and 2nd dimensions of data must have the same sizes')
ntim  = size(dat,3);
assert(ntim==length(opt.timelabel),'opt.timelabel must have the same length as the data');
if ~isempty(opt.plot.AlphaData)
  if isvector(opt.plot.AlphaData)
    assert(size(opt.plot.AlphaData,1)==nchan,'Mask data must consist of the same channels as the data.');
    maskdat = repmat(opt.plot.AlphaData,[1,size(dat,[2:ndims(dat)])]);
  elseif ~ismatrix(opt.plot.AlphaData)
    assert(isequal(size(dat),size(opt.plot.AlphaData)),'Mask data must be the same size as the data.');
    maskdat = opt.plot.AlphaData;
  elseif size(opt.plot.AlphaData,1)==size(opt.plot.AlphaData,2)
    assert(size(opt.plot.AlphaData,1)==nchan,'Mask data must consist of the same channels as the data.');
    maskdat = repmat(opt.plot.AlphaData,[1,1,size(dat,[3:ndims(dat)])]);
  else
    assert(isequal(size(dat,[1,3]),size(opt.plot.AlphaData)),'Mask data must consist of the same channels and the same time points as the data.');
    maskdat = repmat(opt.plot.AlphaData,[1,1,size(dat,[2,4:ndims(dat)])]);
    maskdat = permute(maskdat,[1,3,2,4:ndims(maskdat)]);
  end
else
  maskdat = ones(size(dat));
end
%-- no seed
if nxchan ==1
    opt.seed             = 1;
    opt.plot.Interaction = false;
    opt.plot.ShowSeed    = false;
end

%-- set seed
if isnumeric(opt.seed)
  if numel(opt.seed)==1 && isint(opt.seed)
    seedchidx = opt.seed;
  else
    error('Only one channel can be specified as seed.');
  end
else
    seedchidx = cellstrfind(opt.channels.name,opt.seed);
end

%%% run plotGridSCCommon
plotGridSCCommon;

%-- Plot figures
if isnumeric(opt.plot.Time)
    itim = round(opt.plot.Time);
else
    itim = find(ismember(opt.timelabel,opt.plot.Time),1);
end
  for ee = 1:length(inx)
    %-- Decide how many subplots are needed
    nRow = opt.plot.nSubPlots(1);
    nCol = opt.plot.nSubPlots(2);
    
    plotElectrodes = gridList(inx{ee});
    
    %%% construct image data & alpha data for NaN
    scidx = reshape(inx{ee},nCol,nRow)';
    sensidx = grididx(scidx);
    imdat = nan(nRow,nCol);
    imdat(~isnan(sensidx))   = dat(sensidx(~isnan(sensidx)),seedchidx,itim);
    imalp = ones(nRow,nCol);
    imalp(~isnan(sensidx))   = maskdat(sensidx(~isnan(sensidx)),seedchidx,itim);
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
    if nxchan ==1
        htl = title(opt.timelabel(itim));
    else
        htl = title(sprintf('%s - %s',opt.channels.name{seedchidx},opt.timelabel{itim}));
    end
    figlist{ee} = hF;
    
    %%% Show seed
    hold(hsc,'on');

    %-- Show gray patch & cross
    seeddat = ones([size(him.CData),3])*.6;
    seedgrididx = find(sensidx==seedchidx,1);
    hsd = imagesc(hsc,seeddat);
%     hsd.CDataMapping='direct';
    hsd.AlphaData=zeros(size(hsd.CData,[1,2]));
    hx = [];
        
    if opt.plot.ShowSeed
        if ~isempty(seedgrididx)
            hsd.AlphaData(seedgrididx)=1;
            
            seedgridpos = [ceil(seedgrididx./nRow),mod(seedgrididx-1,nRow)+1];
            hx    = plot(hsc,seedgridpos(1)+[-.5 .5],seedgridpos(2)+[-.5 .5],'k-','LineWidth',1.5);
            hx(2) = plot(hsc,seedgridpos(1)+[-.5 .5],seedgridpos(2)+[.5 -.5],'k-','LineWidth',1.5);
        end

    end
    
    %%% Set UI
    him.UserData.sensidx   = sensidx;
    him.UserData.seedchidx = seedchidx;
    him.UserData.ShowSeed  = opt.plot.ShowSeed;
    him.UserData.timelabel = opt.timelabel;
    him.UserData.channame  = opt.channels.name;
    him.UserData.htl       = htl;
    him.UserData.hsd       = hsd;
    him.UserData.hx        = hx;
    
    him.UserData.dat = dat;
    him.UserData.maskdat = maskdat;
    
    %-- for Update time
    if ntim > 1
        sldpos = hsc.Position * [1,0,0,0;0,-0.25,0,0;0,0,1,0;0,1,0,0]';
        hsl = uicontrol('Style','slider','Units','normalized','Position',sldpos,...
            'value',itim, 'min',1, 'max',ntim);
        hsl.Callback = @(es,ed) UpdateTime(him,hsl);
        hsl.SliderStep = [1/ntim min(5/ntim,1)];

        %-- for Update time by Keybroad
        hF.KeyPressFcn = @(es,ed) UpdateTimeKey(ed,him,hsl);
    else
       hsl = [];
       hsl.Value = 1;
    end
    
    %-- for Update seed
    if opt.plot.Interaction
        hbg.ButtonDownFcn = @(es,ed) UpdateSeed(ed,him,hsl);
    end
    
  end
end

%% UI Functions
function UpdateSC(him, hsl)
    itim      = hsl.Value;
    [nRow, nCol] = size(him.CData);
    sensidx   = him.UserData.sensidx;
    seedchidx = him.UserData.seedchidx;
    timelabel = him.UserData.timelabel;
    channame  = him.UserData.channame;
    
    dat = him.UserData.dat;
    maskdat = him.UserData.maskdat;
    nxchan  = size(dat,2);

    imdat = nan(nRow,nCol);
    imdat(~isnan(sensidx))   = dat(sensidx(~isnan(sensidx)),seedchidx,itim);
    imalp = ones(nRow,nCol);
    imalp(~isnan(sensidx))   = maskdat(sensidx(~isnan(sensidx)),seedchidx,itim);
    imalp(isnan(imdat))   = 0;
    
    him.CData       = imdat;
    him.AlphaData   = imalp;
    
    %-- update title
    htl         = him.UserData.htl;
    htl.String  = timelabel(itim);
    if nxchan ==1
        htl.String = timelabel(itim);
    else
        htl.String = sprintf('%s - %s',channame{seedchidx},timelabel{itim});
    end
end

function UpdateTime(him, hsl)
    hsl.Value   = round(hsl.Value);
    
    %-- Update plot
    UpdateSC(him, hsl);
end

function UpdateTimeKey(ed, him, hsl)
    itim        = hsl.Value;
    getKey      = ed.Key;
    if ismember({'rightarrow'},getKey) && itim<hsl.Max
        shifttim    = 1;
    elseif  ismember({'leftarrow'},getKey) && itim>hsl.Min
        shifttim    = -1;
    else
        shifttim    = 0;
    end
    
    hsl.Value   = round(hsl.Value)+shifttim;
    
    %-- Update plot
    UpdateSC(him, hsl);
end

function UpdateSeed(ed, him, hsl)
    hsd         = him.UserData.hsd;
    hx          = him.UserData.hx;
    ShowSeed    = him.UserData.ShowSeed;
    clickpos    = round(ed.IntersectionPoint(1:2));
    
    %-- Update plot
    sensidx     = him.UserData.sensidx;
    seedchidx   = sensidx(clickpos(2),clickpos(1));
    if ~isnan(seedchidx)
        him.UserData.seedchidx = seedchidx;
        UpdateSC(him, hsl);
        
        %-- show seed pos
        seedgrididx = find(sensidx==seedchidx,1);
        if ShowSeed && ~isempty(seedgrididx)
            hsd.AlphaData(:)=0;
            hsd.AlphaData(seedgrididx)=1;
            if isempty(hx)
                hsc = hsd.Parent;
                hx    = plot(hsc,clickpos(1)+[-.5 .5],clickpos(2)+[-.5 .5],'k-','LineWidth',1.5);
                hx(2) = plot(hsc,clickpos(1)+[-.5 .5],clickpos(2)+[.5 -.5],'k-','LineWidth',1.5);
                him.UserData.hx  = hx;
            else
                hx(1).XData = clickpos(1)+[-.5 .5];
                hx(1).YData = clickpos(2)+[-.5 .5];
                hx(2).XData = clickpos(1)+[-.5 .5];
                hx(2).YData = clickpos(2)+[.5 -.5];
            end
        end
    end
    
end

