function figlist = ecog_plotGridBarwitherr(dat, errlo, errup, whichGrid, opt)

% p = ecog_plotGridBarwitherr(data, lower_error, upper_error, whichGrid, options)
% Plot bars with error bars in the entire grid.
% 
% dat, errlo, errup = datasets consisting of channels x samples matrix
% whichGrid         = 'G', 'GA' or 'GB'
% 
% options mush include 'channels'
% options should include 'timelable','viselec'
%   options.plot.tickCol = 'strings': colors of tick labels
%                        or 
%   options.plot.tickCol = {'strings', 'expression'} : 
%      	You can specify multiple tick colors with expression using @function such as
%           {'red', @(x,y,z) min(x)>=0;  'blue', @(x,y,z) max(x)<=0},
%       where expression takes three arguemnts of "dat(channel,:)",
%       "errlo(channel,:)" and "errup(channel,:)".
%       If multiple expressions are satifsied simultaneously, the first color is applied.
%       If no expression is satifsied, the tick color is 'black'.

% Dependency: plotGridBarCommon, barwitherr, SetDefault, cellstrfind

% 20190726 Yuasa
% 20190813 Yuasa: fix ylim computation
% 20220222 Yuasa: use plotGridBarCommon

%% Set options
narginchk(4,inf);
SetDefault('opt.plot.nSubPlots',[]);
SetDefault('opt.plot.addEccToTitle','no');
SetDefault('opt.plot.fontSize',12);
SetDefault('opt.plot.YLim',[]);
SetDefault('opt.plot.FigName','');
SetDefault('opt.plot.RotGrid',false);
SetDefault('opt.plot.colors',[0:(size(dat,2)-1)]);
SetDefault('opt.plot.tickCol',[]);

% whichHDgrid     = upper(whichHDgrid);
%%
%-- plot color properties
if isvector(dat),   ngroups = 1;
else,               ngroups = size(dat,2);
end
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

%-- set values for error bar
if isempty(errlo),  errlo   = nan(size(dat));   end
if isempty(errup),  errup   = nan(size(dat));   end
iserrlo = ~all(isnan(errlo(:)));
iserrup = ~all(isnan(errup(:)));
nanlo = isnan(errlo);    errlo(nanlo) = 0;
nanup = isnan(errup);    errup(nanup) = 0;

%%% run plotGridBarCommon
tickColChk = 'opt.plot.tickCol{itcl,2}(dat(ichs(1),:),errlo(ichs(1),:),errup(ichs(1),:))';
plotGridBarCommon;

%-- Set ylim
ichs = cellstrfind(opt.channels.name,gridList,1);
if isempty(opt.plot.YLim)
    if iserrlo,    ymin = prctile(reshape(dat(ichs,:)-errlo(ichs,:),1,[]),5);
    else,          ymin = prctile(reshape(dat(ichs,:),1,[]),0.3);     end
    if iserrup,    ymax = prctile(reshape(dat(ichs,:)+errup(ichs,:),1,[]),95);
    else,          ymax = prctile(reshape(dat(ichs,:),1,[]),99.7);     end
    if ~iserrlo,   yminex = (ymin - (ymax-ymin)*0.1);  ymin = yminex .* (ymin<0 | yminex>0);
    end
    if ~iserrup,   ymaxex = (ymax + (ymax-ymin)*0.1);  ymax = ymaxex .* (ymax>0 | ymaxex<0);
    end
    if ymax-ymin>2
        ymin = floor(ymin); ymax = ceil(ymax);
    elseif ymax-ymin>1
        ymin = floor(ymin*2)/2; ymax = ceil(ymax*2)/2;
    else
        ymin = floor(ymin*10)/10; ymax = ceil(ymax*10)/10;
    end
else
    ymin = opt.plot.YLim(1);
    ymax = opt.plot.YLim(2);
end

%-- Legend
flglgnd = isfield(opt,'timelabel');

%-- Plot figures
errlo(nanlo) = nan;    errup(nanup) = nan;
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
      plterrlo = nan(length(plotElectrodes), size(errlo,2));
      plterrup = nan(length(plotElectrodes), size(errup,2));
      pltdat = nan(length(plotElectrodes), size(dat,2));
        plterrlo(grididx,:) = errlo(ichs,:);
        plterrup(grididx,:) = errup(ichs,:);
        pltdat(grididx,:)   = dat(ichs,:);
      figure(hF);  hA = subplot(nRow,1,ichidx);
      bars = barwitherr(cat(3,plterrlo,plterrup), pltdat);
      set(hA,'XTickLabel',nameElectrodes,'XTick',1:numel(plotElectrodes),...
              'XAxisLocation','top','Ylim',[ymin ymax],'FontSize', opt.plot.fontSize,...
              'Position',get(hA,'Position')+[-0.08 0 0.14 0]);
      for ibar=1:length(bars),  bars(ibar).FaceColor = plotCol(ibar,:);  end
      if flglgnd && (ee ==1 || ee ~= length(inx)) && ichidx == 1
          nlegcol = 1 + (ngroups>7);
          legend(opt.timelabel,'Location','best','NumColumns',nlegcol);
      elseif flglgnd && ee ~= 1 && ee == length(inx) && ichidx == nRow
          nlegcol = 1 + (ngroups>7);
          legend(opt.timelabel,'Location','best','NumColumns',nlegcol);
      end
    end
    figlist{ee} = hF;
end

