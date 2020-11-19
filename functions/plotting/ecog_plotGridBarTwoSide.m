function figlist = ecog_plotGridBarTwoSide(dat1, dat2, whichHDgrid, opt, errlo1, errup1, errlo2, errup2)

% p = ecog_plotGridBarTwoSide(data1, data2, whichHDgrid, options, [lower_error1, upper_error1, lower_error2, upper_error2])
% Plot the entire HD grid
% 
% options mush include 'channels'
%   options.plot.tickCol = 'strings': colors of tick labels
%                        {'strings', @function} : change color based on the function 

% Dependency: SetDefault, cellstrfind

% 20200321 Yuasa
% 20200428 Yuasa: enable errorbar

%% Set options
narginchk(4,inf);

haserror = nargin > 4;
SetDefault('errlo1',nan(size(dat1)));
SetDefault('errup1',nan(size(dat1)));
SetDefault('errlo2',nan(size(dat2)));
SetDefault('errup2',nan(size(dat2)));

SetDefault('opt.plot.nSubPlots',[]);
SetDefault('opt.plot.addEccToTitle','no');
SetDefault('opt.plot.fontSize',12);
SetDefault('opt.plot.FigName','');
SetDefault('opt.plot.RotGrid',false);
SetDefault('opt.plot.colors',[1,2]);
SetDefault('opt.plot.tickCol',[]);

SetDefault('opt.plot.YLim',[],'cell');
SetDefault('opt.plot.axis','xy','cell');

if length(opt.plot.YLim)==1
    opt.plot.YLim{2} = [];
end
if length(opt.plot.axis)==1
    opt.plot.axis = repmat(opt.plot.axis,1,2);
end

% whichHDgrid     = upper(whichHDgrid);
%%
%-- plot color properties
ngroups = 2;
if isnumeric(opt.plot.colors)
    if (size(opt.plot.colors,2)~=3 || any(opt.plot.colors(:)>1))
        %%% specify as integer
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
        %%% specify as rgb
        plotCol = opt.plot.colors;
    end
elseif ischar(opt.plot.colors)
    %%% specify as 'r' 'g' 'b'
    plotCol = reshape(opt.plot.colors,[],1);
end 
plotCol = repmat(plotCol,ceil(ngroups./size(opt.plot.colors,1)),1);
plotCol = plotCol(1:ngroups,:);

%-- correct elecnames (set '%03d')
[eleccat, elecnum] = strtok(opt.channels.name,int2str(0:9));
for el = 1:length(opt.channels.name)
    opt.channels.name{el} = sprintf('%s%03d',eleccat{el},str2double(elecnum{el}));
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

%-- Create a list of electrode names for the overall grid layout
gridList = [];
areaList = [];
isTableCol = @(t, thisCol) ismember(thisCol, t.Properties.VariableNames);
for ee = 1:FullGRID
    gridList{ee} = sprintf('%s%03d',whichHDgrid,ee);
    
    %%% get Visual Area Label
    areaList{ee} = sprintf('%s',gridList{ee});
    igrid = find(strcmp(opt.channels.name,gridList{ee}),1);
    if ~isempty(igrid)
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
        
      if nCol < 12 % avoid crowded tick labels
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
      end
    end
    TickCol = 'black';
    if ~isempty(opt.plot.tickCol)
        if numel(opt.plot.tickCol)==1, TickCol = char(opt.plot.tickCol);
        else
          assert(size(opt.plot.tickCol,2)==2, 'specs.plot.tickCol is invalid');
          ichs = cellstrfind(opt.channels.name,gridList{ee});
          if ~isempty(ichs)
            for itcl=size(opt.plot.tickCol,1):-1:1
              if opt.plot.tickCol{itcl,2}(dat1(ichs(1),:),dat1(ichs(1),:))
                TickCol = char(opt.plot.tickCol{itcl,1});
              end
            end
          end
        end
    end
    areaList{ee} = sprintf('{\\color{%s} \\bf %s}',TickCol,areaList{ee});
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

if isempty(opt.plot.YLim{1})
    y1min = nanmin(reshape(dat1(ichs,:),1,[]));
    y1max = nanmax(reshape(dat1(ichs,:),1,[]));
    
    ystep = (y1max - y1min)./5;
    ystep = round(ystep,-floor(log10(ystep)));
    
    y1min = floor(y1min./ystep).*ystep;
    y1max = ceil(y1max./ystep).*ystep;
else
    y1min = opt.plot.YLim{1}(1);
    y1max = opt.plot.YLim{1}(2);
end
yyl = [y1min y1max];
if ismember(opt.plot.axis(1),'ij')
    yyl = fliplr(-yyl);
end
if isempty(opt.plot.YLim{2})
    y2min = nanmin(reshape(dat2(ichs,:),1,[]));
    y2max = nanmax(reshape(dat2(ichs,:),1,[]));
    
    ystep = (y2max - y2min)./5;
    ystep = round(ystep,-floor(log10(ystep)));
    
    yyr = [floor(y2min./ystep), ceil(y2max./ystep)].*ystep;
    if ismember(opt.plot.axis(2),'ij')
        yyr = fliplr(-yyr);
    end
    
    lrratio = yyr./yyl;
    if lrratio(1)>lrratio(2)
        yyr(2) = yyl(2)./yyl(1).*yyr(1);
    else
        yyr(1) = yyl(1)./yyl(2).*yyr(2);
    end
    if ismember(opt.plot.axis(2),'ij')
        yyr = fliplr(-yyr);
    end
    
    y2min = yyr(1);
    y2max = yyr(2);
else
    y2min = opt.plot.YLim{2}(1);
    y2max = opt.plot.YLim{2}(2);
end

%-- Legend
flglgnd = isfield(opt.plot,'legend');

%-- Plot figures
for ee = 1:length(inx)
    %-- Decide how many subplots are needed
    nRow = opt.plot.nSubPlots(1);
    nCol = opt.plot.nSubPlots(2);

    %%% Bar plot
    hF = figure('Name', opt.plot.FigName); 
    set(hF, 'Position', [150 100 2000 1250]);
    for ichidx=1:nRow
      plotElectrodes = gridList(inx{ee}((ichidx-1)*nCol + (1:nCol)));
      nameElectrodes = areaList(inx{ee}((ichidx-1)*nCol + (1:nCol)));
      ichs = cellstrfind(opt.channels.name,plotElectrodes,0);
      grididx = cellstrfind(plotElectrodes,opt.channels.name(ichs),1);
      pltdat1 = nan(length(plotElectrodes), size(dat1,2));
      pltdat2 = nan(length(plotElectrodes), size(dat2,2));
        pltdat1(grididx,:)   = dat1(ichs,:);
        pltdat2(grididx,:)   = dat2(ichs,:);
      set(0,'CurrentFigure',hF);  hA = subplot(nRow,1,ichidx);
      
      yyaxis left
      bars = bar((1:nCol)-0.2,pltdat1,0.4);       % bar plot
      set(hA,'XTickLabel',nameElectrodes,'XTick',1:numel(plotElectrodes),...
              'XAxisLocation','top','Ylim',[y1min y1max],'FontSize', opt.plot.fontSize,...
              'Position',get(hA,'Position')+[-0.08 0 0.14 0]);
      axis(hA,opt.plot.axis{1});
          
      yyaxis right
      bars(2) = bar((1:nCol)+0.2,pltdat2,0.4);       % bar plot
      set(hA,'Ylim',[y2min y2max]);
      axis(hA,opt.plot.axis{2});
          
      for ibar=1:length(bars),  bars(ibar).FaceColor = plotCol(ibar,:);  end
      if flglgnd && (ee ==1 || ee ~= length(inx)) && ichidx == 1
          nlegcol = 1 + (ngroups>7);
          legend(opt.plot.legend,'Location','best','NumColumns',nlegcol);
      elseif flglgnd && ee ~= 1 && ee == length(inx) && ichidx == nRow
          nlegcol = 1 + (ngroups>7);
          legend(opt.plot.legend,'Location','best','NumColumns',nlegcol);
      end
    end
    figlist{ee} = hF;
end

