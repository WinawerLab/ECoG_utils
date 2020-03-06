function trials_out = ecog_plotGridPRF(whichHDgrid, opt, varargin)

% trials_out = ECOG_PLOTGRIDPRF(whichHDgrid, option, data1)
% trials_out = ECOG_PLOTGRIDPRF(whichHDgrid, option, data1, data2,...)
% ECOG_PLOTGRIDPRF draw pRF in the entire HD grid.
% 
% data is a structure of pRF information, requiring: 
%   data.ecc:    Nx1 array with pRF eccentricity
%   data.ang:    Nx1 array with pRF angle [degree: 0deg is on positive x-axis and increase anti-clockwise]
%   data.rfsize: Nx1 array with pRF size
% each data can be Mx1 cell-array with pRF structures.
% 
% option is a structure of plot information.
%   option.channels     = a table of channel information.
%  	option.viselec      = (option) a structure of visual area information.
%   option.plot         = (option) a structure of plot information with following fields.
%       Type            = 'contour'(default) or 'imagesc'
%       sigma           = scalar | vector (default = 1)
%                         To draw contour lines
%                           - from s.d. to N*s.d., specify as the scalar value N.
%                           - at specific s.d.s, specify as a vector of coefficients of s.d.s.
%                             Specify as [1:N] is the same as [N].
%                           - at a specific K*s.d., specify as a two-element row vector [K K].
%       isavg           = true or flase (default), if true plot averaged prf within each data.
%       resolution      =
%       cfactor         =
%       XLim            =
%       YLim            =
%       colors          =
% 
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
% Followings are unexplained options for option.plot: 
%       representative
%       rep_options

% option should include 'channels'

% Dependency: SetDefault, cellstrfind

% 20191106 Yuasa
% 20200206 Yuasa: change arguments structures (not compatible to prior version)

%% Set options
narginchk(3,inf);
assert(isfield(opt,'channels')&&~isempty(opt.channels),'options.channels is required');

SetDefault('opt.plot.Type', 'contour', 0);
SetDefault('opt.plot.isavg', false, 0);
SetDefault('opt.plot.addEccToTitle', 'no', 0);
SetDefault('opt.plot.fontSize', 12, 0);
SetDefault('opt.plot.sigma', [1], 1);
SetDefault('opt.plot.resolution', 100, 1);
SetDefault('opt.plot.cfactor', 1, 1);
SetDefault('opt.plot.nSubPlots', [], 1);
SetDefault('opt.plot.XLim', [], 1);
SetDefault('opt.plot.YLim', [], 1);

SetDefault('opt.plot.FigName', '', 0);
SetDefault('opt.plot.RotGrid', false, 0);
SetDefault('opt.plot.colors', [], 1);
SetDefault('opt.plot.labelCol', [], 1);
SetDefault('opt.plot.representative', [], 1);
SetDefault('opt.plot.legend', {}, 1);
SetDefault('opt.plot.options', {}, 1);
SetDefault('opt.plot.rep_options', {}, 1);

whichHDgrid     = upper(whichHDgrid);
%%
%-- data correction (align for multi group)
ngroups = nargin-2;
groups  = [];
for ii=1:ngroups
    if iscell(varargin{ii}),  varargin{ii} = reshape(varargin{ii},[],1);
    else,                     varargin{ii} = varargin(ii);
    end
    groups  = [groups; ii*ones(size(varargin{ii}))];
end
dat      = cat(1,varargin{:});
ndats    = numel(dat);

%-- representative data
if ~iscell(opt.plot.representative)
    opt.plot.representative = repmat({reshape(opt.plot.representative,[],1)},ngroups,1);
end
repidx = false(0);
for ii=1:ngroups
    if islogical(opt.plot.representative{ii})
        tmpidx = reshape(opt.plot.representative{ii},[],1);
    else
        tmpidx = false(numel(find(groups==ii)),1);
        tmpidx(opt.plot.representative{ii}) = true;
    end
    repidx = [repidx; tmpidx];
end

%-- plot options
if numel(opt.plot.options)~=ngroups || ~any(cellfun(@iscell,opt.plot.options))
    opt.plot.options = repmat({opt.plot.options},ngroups);
end
if numel(opt.plot.rep_options)~=ngroups || ~any(cellfun(@iscell,opt.plot.rep_options))
    opt.plot.rep_options = repmat({opt.plot.rep_options},ngroups);
end

%-- set s.d. to plot
if isscalar(opt.plot.sigma)
    if opt.plot.sigma>=2,   opt.plot.sigma = 1:opt.plot.sigma;
    else,                   opt.plot.sigma = [opt.plot.sigma opt.plot.sigma];
    end
end

%-- correct elecnames (set '%03d')
[eleccat, elecnum] = strtok(opt.channels.name,int2str(0:9));
for el = 1:length(opt.channels.name)
    opt.channels.name{el} = sprintf('%s%03d',eleccat{el},str2double(elecnum{el}));
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
if ngroups>1
    plotCol = repmat(plotCol,ceil(ngroups./size(plotCol,1)),1);
    plotCol = plotCol(1:ngroups,:);
    reptCol = plotCol;
else
    plotCol = repmat(plotCol,ceil(ndats./size(plotCol,1)),1);
    plotCol = plotCol(1:ndats,:);
    reptCol = 'r';
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
                    if isTableCol(opt.channels,'bensoneccen')
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

%-- get gaussian properties of pRF & group updates
grididx = cellstrfind(opt.channels.name,gridList,1);
PRFs = cell(size(dat));
for ii = 1:ndats
    PRFs{ii} = [dat{ii}.ecc(grididx,1) .* cos(dat{ii}.ang(grididx,1)/180*pi), ...
                dat{ii}.ecc(grididx,1) .* sin(dat{ii}.ang(grididx,1)/180*pi), ...
                dat{ii}.rfsize(grididx,1)];
end
PRFs = cat(3,PRFs{:});
PRFr = PRFs(:,:,repidx);
PRFs = PRFs(:,:,~repidx);
groupr = groups(repidx);
groups = groups(~repidx);
%%-- PRFa is determined to estimate proper axis size
if opt.plot.isavg
    PRFa = [];
    for ii=1:ngroups
        PRFa = cat(3,nanmean(PRFs(:,:,groups==ii), 3));
    end
else
    PRFa = PRFs;
    if ngroups > 1
        plotCol = plotCol(groups,:);
    end
end
PRFa = cat(3,PRFa,PRFr,PRFr);

%-- legend
lgnds = reshape(opt.plot.legend,1,[]);
if ~iscell(lgnds), lgnds = {lgnds}; end
flglgnd = ~isempty(lgnds);
if ngroups>1
    lgndr = lgnds([1:length(lgnds)]>ngroups);
    lgnds = lgnds([1:length(lgnds)]<=ngroups);
else
    lgndr = lgnds(find(repidx)<=length(lgnds));
    lgnds = lgnds(find(~repidx)<=length(lgnds));
end

%-- set ylim
if ~isempty(opt.plot.XLim) && ~isempty(opt.plot.YLim)
        xmin = opt.plot.XLim(1) ./opt.plot.cfactor;
        xmax = opt.plot.XLim(2) ./opt.plot.cfactor;
        ymin = opt.plot.YLim(1) ./opt.plot.cfactor;
        ymax = opt.plot.YLim(2) ./opt.plot.cfactor;
else
    xmin = prctile(prctile(PRFa(:,1,:)-max(opt.plot.sigma).*PRFa(:,3,:),5,3),5,1);
    xmax = prctile(prctile(PRFa(:,1,:)+max(opt.plot.sigma).*PRFa(:,3,:),95,3),95,1);
    ymin = prctile(prctile(PRFa(:,2,:)-max(opt.plot.sigma).*PRFa(:,3,:),5,3),5,1);
    ymax = prctile(prctile(PRFa(:,2,:)+max(opt.plot.sigma).*PRFa(:,3,:),95,3),95,1);
    if ~isempty(opt.plot.YLim)
        xmin = opt.plot.XLim(1) ./opt.plot.cfactor;
        xmax = opt.plot.XLim(2) ./opt.plot.cfactor;
        axlen = xmax - xmin;
        axcnt = mean([ymin,ymax]);
        ymin = round(axcnt,-round(log10(axlen)-2)) - axlen/2;
        ymax = round(axcnt,-round(log10(axlen)-2)) + axlen/2;
    elseif ~isempty(opt.plot.XLim)
        ymin = opt.plot.YLim(1) ./opt.plot.cfactor;
        ymax = opt.plot.YLim(2) ./opt.plot.cfactor;
        axlen = ymax - ymin;
        axcnt = mean([xmin,xmax]);
        xmin = round(axcnt,-round(log10(axlen)-2)) - axlen/2;
        xmax = round(axcnt,-round(log10(axlen)-2)) + axlen/2;
    else
        axlen = max(xmax-xmin,ymax-ymin);
        axlen = round(axlen,-round(log10(axlen)-1));
        axcnt = mean([xmin,xmax]);
        xmin = round(axcnt,-round(log10(axlen)-2)) - axlen/2;
        xmax = round(axcnt,-round(log10(axlen)-2)) + axlen/2;
        axcnt = mean([ymin,ymax]);
        ymin = round(axcnt,-round(log10(axlen)-2)) - axlen/2;
        ymax = round(axcnt,-round(log10(axlen)-2)) + axlen/2;
    end
end

%-- main
axmaxmin = max(xmax,ymax)-min(xmin,ymin);
axlen    = max(xmax-xmin,ymax-ymin);
axspd    = axmaxmin./axlen;
axplt    = linspace(min(xmin,ymin),max(xmax,ymax),ceil(opt.plot.resolution * axspd));
[xi, yi]=meshgrid(axplt,axplt);
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
          %-- prepare gaussian image
          PRFexp = (bsxfun(@minus,xi,PRFs(el,1,:)).^2 + bsxfun(@minus,yi,PRFs(el,2,:)).^2)./(2.*PRFs(el,3,:).^2);
          PRFamp = 1;
%               PRFamp = 1 ./ (PRFs(el,3,ii) * sqrt(2*pi) .* opt.plot.cfactor);
          imPRFs = PRFamp  .* exp(-PRFexp);
          if opt.plot.isavg
              imPRFa = [];
              for ii=1:ngroups,    imPRFa = cat(3,imPRFa,nanmean(imPRFs(:,:,groups==ii),3));    end
              imPRFs = imPRFa;
          end
          PRFexp = (bsxfun(@minus,xi,PRFr(el,1,:)).^2 + bsxfun(@minus,yi,PRFr(el,2,:)).^2)./(2.*PRFr(el,3,:).^2);
          PRFamp = 1;
%               PRFamp = 1 ./ (PRFr(el,3,ii) * sqrt(2*pi) .* opt.plot.cfactor);
          imPRFr = PRFamp  .* exp(-PRFexp);
                    
          %-- plot       hA = subplot(nRow,nCol,ichidx);
          set(0,'CurrentFigure',hF);
          ax_sp = [0.03 0.03]; ax_siz = (1-ax_sp*2) ./[nCol nRow];
          ax_pos = [ax_sp(1) + ax_siz(1)*(mod(ichidx-1,nCol)),...
                    ax_sp(2) + ax_siz(2)*(nRow-ceil(ichidx./nCol)),...
                    ax_siz(1)*0.70, ax_siz(2)*0.70];
          hA = subplot('Position',ax_pos);    daspect([1 1 1]);
              outerpos = hA.OuterPosition;
              ax_ti = hA.TightInset;
              ax_left = outerpos(1) + ax_ti(1);
              ax_bottom = outerpos(2) + ax_ti(2);
              ax_width = outerpos(3) - ax_ti(1) - ax_ti(3);
              ax_height = outerpos(4) - ax_ti(2) - ax_ti(4);
              hA.Position = [ax_left ax_bottom ax_width ax_height];
          hold on
          plot([min(axplt) max(axplt)],[0 0],'k:');
          plot([0 0],[min(axplt) max(axplt)],'k:');
          hplts = [];
          switch opt.plot.Type
              case {'imagesc'}
                  imagesc(axplt.*opt.plot.cfactor,axplt.*opt.plot.cfactor,zeros(length(axplt)),...
                              opt.plot.options{groups(ii)}{:},[0 exp(-0.5.^2./2)]);
                  axis xy;
          end
          for ii = 1:size(imPRFs,3)
              switch opt.plot.Type
                  case {'contour'}
                      [~,hplts(ii)] = contour(axplt.*opt.plot.cfactor,axplt.*opt.plot.cfactor,imPRFs(:,:,ii),...
                                  exp(-opt.plot.sigma.^2./2), 'LineColor', plotCol(ii,:),'LineWidth',0.5,...
                                  opt.plot.options{groups(ii)}{:});
                  case {'imagesc'}
                      hplts(ii) = imagesc(axplt.*opt.plot.cfactor,axplt.*opt.plot.cfactor,imPRFs(:,:,ii),...
                                  'AlphaData',ones(size(imPRFs,[1 2]))./(size(imPRFs,3)+2*size(imPRFr,3)),...
                                  opt.plot.options{groups(ii)}{:},[0 exp(-0.5.^2./2)]);
                      axis xy;
              end
          end
          hpltr = [];
          for ii = 1:size(imPRFr,3)
              switch opt.plot.Type
                  case {'contour'}
                      [~,hpltr(ii)] = contour(axplt.*opt.plot.cfactor,axplt.*opt.plot.cfactor,imPRFr(:,:,ii),...
                                  exp(-opt.plot.sigma.^2./2),'LineColor', reptCol(groupr(ii),:),'LineWidth',2.0,...
                                  opt.plot.options{groupr(ii)}{:});
                  case {'imagesc'}
                      hpltr(ii) = imagesc(axplt.*opt.plot.cfactor,axplt.*opt.plot.cfactor,imPRFr(:,:,ii),...
                                  'AlphaData',ones(size(imPRFs,[1 2]))./(size(imPRFs,3)+2*size(imPRFr,3))*2,...
                                  opt.plot.rep_options{groupr(ii)}{:},[0 exp(-0.5.^2./2)]);
                      axis xy;
              end
          end
          hold off
          set(hA,'XLim',[xmin xmax].*opt.plot.cfactor,'YLim',[ymin ymax].*opt.plot.cfactor,...
              'FontSize', opt.plot.fontSize);
          title(nameElectrodes);
          if flglgnd
              lgndsidx = ~cellfun(@isempty,lgnds);      lgndridx = ~cellfun(@isempty,lgndr);
              hlg = legend([hplts(lgndsidx) hpltr(lgndridx)],[lgnds(lgndsidx) lgndr(lgndridx)],...
                  'Location','southwest','FontSize',opt.plot.fontSize);
              hlg.Position = hlg.Position .* [0 1 1 1] + [hA.OuterPosition(1)-hlg.Position(3)-0.01 0 0 0];
              flglgnd = false;
          end
%       if (ee ==1 || ee ~= length(inx)) && ichidx == 1
%           nlegcol = 1 + (ngroups>7);
%           legend(opt.timelabel,'Location','best','NumColumns',nlegcol);
%       elseif ee ~= 1 && ee == length(inx) && ichidx == nRow
%           nlegcol = 1 + (ngroups>7);
%           legend(opt.timelabel,'Location','best','NumColumns',nlegcol);
%       end
          drawnow limitrate nocallbacks;
      end
    end
    trials_out{ee} = hF;
end
