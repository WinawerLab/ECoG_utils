function figlist = ecog_plotGridBarTwoSide(dat1, dat2, whichGrid, opt, errlo1, errup1, errlo2, errup2)

% p = ecog_plotGridBarTwoSide(data1, data2, whichGrid, options, [lower_error1, upper_error1, lower_error2, upper_error2])
% Plot bars for two datasets in the entire grid.
% 
% Input: 
% dat1, dat2    = datasets consisting of channels x samples matrix
% whichGrid     = 'G', 'GA' or 'GB'
% 
% options mush include 'channels'
%   options.plot.tickCol = 'strings': colors of tick labels
%                        or 
%   options.plot.tickCol = {'strings', 'expression'} : 
%      	You can specify multiple tick colors with expression using @function such as
%           {'red', @(x,y) all(x>=y);  'blue', @(x,y) all(x<=y)},
%       where expression takes two arguemnts of "dat1(channel,:)" and "dat2(channel,:)".
%       If multiple expressions are satifsied simultaneously, the first color is applied.
%       If no expression is satifsied, the tick color is 'black'.

% Dependency: plotGridBarCommon, barwitherr, SetDefault, cellstrfind

% 20200321 Yuasa
% 20200428 Yuasa: enable errorbar
% 20220222 Yuasa: use plotGridBarCommon

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

%%% run plotGridBarCommon
tickColChk = 'opt.plot.tickCol{itcl,2}(dat1(ichs(1),:),dat1(ichs(1),:))';
plotGridBarCommon;
 
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
      if haserror
          plterrlo1 = nan(length(plotElectrodes), size(errlo1,2));
          plterrlo2 = nan(length(plotElectrodes), size(errlo2,2));
          plterrup1 = nan(length(plotElectrodes), size(errup1,2));
          plterrup2 = nan(length(plotElectrodes), size(errup2,2));
            plterrlo1(grididx,:)   = errlo1(ichs,:);
            plterrlo2(grididx,:)   = errlo2(ichs,:);
            plterrup1(grididx,:)   = errup1(ichs,:);
            plterrup2(grididx,:)   = errup2(ichs,:);
      end
      set(0,'CurrentFigure',hF);  hA = subplot(nRow,1,ichidx);
      
      yyaxis left
      %-- bar plot
      if haserror
          bars = barwitherr(cat(3,plterrlo1,plterrup1), (1:nCol)-0.2,pltdat1,0.4);
      else
          bars = bar((1:nCol)-0.2,pltdat1,0.4);
      end
      set(hA,'XTickLabel',nameElectrodes,'XTick',1:numel(plotElectrodes),...
              'XAxisLocation','top','Ylim',[y1min y1max],'FontSize', opt.plot.fontSize,...
              'Position',get(hA,'Position')+[-0.08 0 0.14 0]);
      axis(hA,opt.plot.axis{1});
          
      yyaxis right
      %-- bar plot
      if haserror
          bars(2) = barwitherr(cat(3,plterrlo2,plterrup2), (1:nCol)+0.2,pltdat2,0.4);
      else
          bars(2) = bar((1:nCol)+0.2,pltdat2,0.4);
      end
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

