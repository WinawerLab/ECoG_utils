function trials_out = ecog_plotGridBarwitherr(dat, errlo, errup, whichHDgrid, opt)

% trials_out = ecog_plotGridBarwitherr(data, error-upper, error-lower, whichHDgrid, options)
% Plot the entire HD grid
% 
% options should include 'timelable','channels','viselec'
%   options.plot.tickCol = 'strings': colors of tick labels
%                        {'strings', @function} : change color based on the function 

% Dependency: SetDefault, cellstrfind

% 20190726 Yuasa
% 20190813 Yuasa: fix ylim computation

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

whichHDgrid     = upper(whichHDgrid);
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
              if opt.plot.tickCol{itcl,2}(dat(ichs(1),:),errlo(ichs(1),:),errup(ichs(1),:))
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

%-- Plot figures
errlo(nanlo) = nan;    errup(nanup) = nan;
for ee = 1:length(inx)
    %%% Decide how many subplots are needed
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
      if (ee ==1 || ee ~= length(inx)) && ichidx == 1
          nlegcol = 1 + (ngroups>7);
          legend(opt.timelabel,'Location','best','NumColumns',nlegcol);
      elseif ee ~= 1 && ee == length(inx) && ichidx == nRow
          nlegcol = 1 + (ngroups>7);
          legend(opt.timelabel,'Location','best','NumColumns',nlegcol);
      end
    end
    trials_out{ee} = hF;
end

