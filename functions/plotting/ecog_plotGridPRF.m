function figlist = ecog_plotGridPRF(whichHDgrid, opt, varargin)

% p = ECOG_PLOTGRIDPRF(whichHDgrid, option, result1)
% p = ECOG_PLOTGRIDPRF(whichHDgrid, option, result1, result2,...)
% p = ECOG_PLOTGRIDPRF(whichElectrodes, option, result1)
% ECOG_PLOTGRIDPRF draw pRF in the entire HD grid.
% 
% <whichHDgrid> can be 'GA' or 'GB'.
% Otherwise it wll be interpreted as electrodes names (some options will be ignored).
% 
% <result> is a structure of pRF information, requiring: 
%   result.ecc:    Nx1 array with pRF eccentricity
%   result.ang:    Nx1 array with pRF angle [degree: 0deg is on positive x-axis and increase anti-clockwise]
%   result.rfsize: Nx1 array with pRF size
% each result can be Mx1 cell-array with pRF structures.
% 
% <option> is a structure of plot information.
%   option.channels      = a table of channel information.
%                          If result has channel field, it's not necessary.
%  	option.viselec       = a structure of visual area information.
%                          If channels has visual area table or you don't want to show area name,
%                          it's not necessary.
%   option.plot          = a structure of plot information with following fields.
%       Type             = 'contour'(default) or 'imagesc'
%       sigma            = scalar | vector (default = 1)
%                          To draw contour lines
%                            - from s.d. to N*s.d., specify as the scalar value N.
%                            - at specific s.d.s, specify as a vector of coefficients of s.d.s.
%                              Specify as [1:N] is the same as [N].
%                            - at a specific K*s.d., specify as a two-element row vector [K K].
%       isavg            = true or flase (default), if true plot averaged prf within each result.
%       resolution       = stimulus size in pixel (default = 100);
%       pix2deg          = degree/pixel (default = 1)
%       XLim             =
%       YLim             =
%       colors           =
% 
%       addChsToTitle    = 'yes'(default) or 'no'
%       addWangToTitle   = 'yes'(default) or 'no'
%       addBensonToTitle = 'yes'(default) or 'no'
%       addR2ToTitle     = 'yes' or 'no'(default)  % show R2 or xval of result1
%       addEccToTitle    = 'yes' or 'no'(default)  % show eccentricity of result1
%       addSbjToTile     = 'yes' or 'no'(default)  % show subject name
%       nSubPlots        = [m n], plot m-by-n grid electrodes in a figure.
%       RotGrid          = true or false (default), if true, electrodes are
%                          aligned along row instead of column.
%       labelCol         = color or cell-array of colors, color of the titles
%                          in each axis. defualt is black.
%       FigName          = string, name of the figure.
%       fontSize         = scalar, font size in the figure.
%       legend           = string or cell-array of strings, legend labels.
%       options          = any other options for each axis.
%   option.plotbenson    = 'yes' or 'no'(default)  % plot results estimated from benson atlas
%   option.subjet        = string, subject name is helpful to interpret benson prf.
% 
% Followings are unexplained options for option.plot: 
%       representative
%       rep_options

% option should include 'channels'

% Dependency: <analyzePRF>, SetDefault, cellstrfind, arrangeinrect, istablefield

% 20191106 Yuasa
% 20200206 Yuasa: change arguments structures (not compatible to prior version)
% 20200309 Yuasa: enable to plot no HD grid data
% 20200331 Yuasa: add 'plotbenson' option
% 20200925 Yuasa: solve problem using strings in channels
% 20201110 Yuasa: bug fix

%% Set options
narginchk(3,inf);
if isfield(opt,'channels')&&~isempty(opt.channels)
    channels = opt.channels;
elseif isfield(varargin{1},'channels')&&~isempty(varargin{1}.channels)
    channels = varargin{1}.channels;
else
    error('data or option must have ''channels'' field');
end

SetDefault('opt.plotbenson', 'no', 0);

SetDefault('opt.plot.Type', 'contour', 0);
SetDefault('opt.plot.isavg', false, 0);
SetDefault('opt.plot.addChsToTitle', 'yes', 0);
SetDefault('opt.plot.addWangToTitle', 'yes', 0);
SetDefault('opt.plot.addBensonToTitle', 'yes', 0);
SetDefault('opt.plot.addR2ToTitle', 'no', 0);
SetDefault('opt.plot.addEccToTitle', 'no', 0);
SetDefault('opt.plot.addSbjToTitle', 'no', 0);
SetDefault('opt.plot.fontSize', 12, 0);
SetDefault('opt.plot.sigma', [1], 1);
SetDefault('opt.plot.resolution', 100, 1);
SetDefault('opt.plot.pix2deg', 1, 1);
SetDefault('opt.plot.nSubPlots', [], 1);
SetDefault('opt.plot.XLim', [], 1);
SetDefault('opt.plot.YLim', [], 1);

SetDefault('opt.plot.FigName', '', 0);
SetDefault('opt.plot.RotGrid', true, 0);
SetDefault('opt.plot.colors', [], 1);
SetDefault('opt.plot.labelCol', [], 1);
SetDefault('opt.plot.representative', [], 1);
SetDefault('opt.plot.legend', {}, 1,'cell');
SetDefault('opt.plot.options', {}, 1);
SetDefault('opt.plot.rep_options', {}, 1);

SetDefault('whichHDgrid','*',0);

%-- imagesc can not have a legend
switch opt.plot.Type
    case 'imagesc'
        opt.plot.legend = {};
end

%%
%-- Define variables
plotbenson   = false;
if strcmpi(opt.plotbenson,'yes')
    if all(ismember({'bensoneccen','bensonangle','bensonsigma'},channels.Properties.VariableNames))
        plotbenson   = true;
    else
        warning('Benson atlas information is not found');
    end
end
if ~istablefield(channels,'subject_name')
    if isfield(opt,'subject'),        channels.subject_name(:) = opt.subject;
    else
        if isfield(varargin{1},'subject'), channels.subject_name(:) = varargin{1}.subject;
        else,                              channels.subject_name(:) = '';
        end
    end
end
subjects = cellstr(channels.subject_name);

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

%-- Benson reference
if plotbenson
    subj_L = {'beilen','som648','som692','som708','som748'}; % left hemisphere
    subj_R = {'chaam','som661','som674','som718','som726'}; % right hemisphere
    subj_B = {'som723'}; % bilateral
    assert(~isempty(setdiff(subjects,[subj_L,subj_R,subj_B])),...
            'data includes unknown subject');
        
    %-- detect channels in left hemisphere
    islefths = ismember(subjects,subj_L) | ...
                (ismember(subjects,subj_B) & startsWith(channels.name,'L'));
        
    ngroups = ngroups+1;
    groups  = [groups; ngroups];
    ndats   = ndats+1;
    
    dat{ndats}.ecc      = channels.bensoneccen./opt.plot.pix2deg;
    dat{ndats}.ang      = ((-1).^islefths).*channels.bensonangle+90;
    dat{ndats}.rfsize   = channels.bensonsigma./opt.plot.pix2deg;
end

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

%-- plot color properties
ColList = get(groot,'DefaultAxesColorOrder');
if isempty(opt.plot.colors)
    opt.plot.colors = reshape(1:size(ColList,1),[],1);
end
if isnumeric(opt.plot.colors)
    if (size(opt.plot.colors,2)~=3 || any(opt.plot.colors(:)>1))
        %-- specify as integer
        plotCol = ColList(mod(opt.plot.colors(:)-1,size(ColList,1))+1,:);
        plotCol(opt.plot.colors==0,:) = 0.3;   % special color [0.3 0.3 0.3] for '0'
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

%-- check HD grid
if ischar(whichHDgrid) && ismember(whichHDgrid,{'GA','GB'})
    hasHDgrid   = true;
    gridList    = [];
    chanIdx     = [];
else
    hasHDgrid   = false;
    if islogical(whichHDgrid)
      chanIdx     = find(whichHDgrid);
    elseif isnumeric(whichHDgrid) && isequaln(whichHDgrid,round(whichHDgrid))
      chanIdx     = whichHDgrid;
    else    % list of channel names
      chanIdx     = cellstrfind(channels.name, whichHDgrid,0);
    end
    gridList    = channels.name(chanIdx);
    whichHDgrid = 'Channels';
end

%-- correct elecnames (set '%03d')
if hasHDgrid
[eleccat, elecnum] = strtok(channels.name,int2str(0:9));
for el = 1:length(channels.name)
    channels.name{el} = sprintf('%s%03d',eleccat{el},str2double(elecnum{el}));
end
end

%-- set grid parameters
switch whichHDgrid
    case 'GA',  nCol = 8; nRow = 8;
    case 'GB',  nCol = 8; nRow = 16;
    otherwise,  [nCol,nRow] = arrangeinrect(numel(gridList),2,[4,1]);
end
if hasHDgrid && opt.plot.RotGrid, tmp = nCol; nCol = nRow; nRow = tmp; end
FullGRID = nCol * nRow;

if isempty(opt.plot.nSubPlots)
    nFig = ceil(nRow / nCol ./ 1.2);
    opt.plot.nSubPlots = [ceil(nRow./nFig), nCol];
end
plRow  = opt.plot.nSubPlots(1);
nFig = ceil(nRow./plRow);

%-- Check cross-validation
if isfield(dat{1},'xval') && ~all(isnan(dat{1}.xval))
    R2field = 'xval';
elseif isfield(dat{1},'aggregatedtestperformance') && ~all(isnan(dat{1}.aggregatedtestperformance))
    R2field = 'aggregatedtestperformance';
else
    R2field = 'R2';
end

%-- Create a list of electrode names for the overall grid layout
areaList = [];
isTableCol = @(t, thisCol) ismember(thisCol, t.Properties.VariableNames);
for ee = 1:FullGRID
    if hasHDgrid
      gridList{ee} = sprintf('%s%03d',whichHDgrid,ee);
      igrid = find(strcmp(channels.name,gridList{ee}),1);
      if isempty(igrid),    chanIdx(ee)  = nan;
      else,                 chanIdx(ee)  = igrid(1);
      end
    elseif ee > numel(gridList)
      gridList{ee} = 'none';
      chanIdx(ee)  = nan;
    end
    
    %-- get Visual Area Label
    igrid = chanIdx(ee);
    if strcmpi(opt.plot.addR2ToTitle,'yes') && ~isnan(igrid)
        if strcmpi(opt.plot.addChsToTitle,'yes')
            areaList{ee} = sprintf('%s (%.1f%%)',gridList{ee},dat{1}.(R2field)(igrid));
        else
            areaList{ee} = sprintf('(%.1f%%)',gridList{ee},dat{1}.(R2field)(igrid));
        end
    elseif strcmpi(opt.plot.addChsToTitle,'yes')
        areaList{ee} = sprintf('%s',gridList{ee});
    else
        areaList{ee} = '';
    end
    if (strcmpi(opt.plot.addWangToTitle,'yes')||strcmpi(opt.plot.addBensonToTitle,'yes'))&& ...
            (strcmpi(opt.plot.addChsToTitle,'yes')||strcmpi(opt.plot.addR2ToTitle,'yes'))
        areaList{ee} = sprintf('%s\n',areaList{ee});
    end
    if ~isnan(igrid)
        if strcmpi(opt.plot.addSbjToTitle,'yes') && ~isempty(subjects{ee})
            areaList{ee} = sprintf('%s %s',subjects{igrid}, areaList{ee});
        end
        
        VA = [];
        if isTableCol(channels,'wangarea')
           if iscell(channels.wangarea),    VA = channels.wangarea{igrid};
           else,                            VA = channels.wangarea(igrid);
           end
        elseif isfield(opt,'viselec') && isfield(opt.viselec,'wang15_mplbl')
           vigrid = find(strcmp(opt.viselec.wang15_mplbl.elec_labels,gridList{ee}),1);
           if ~isempty(vigrid),  VA = opt.viselec.wang15_mplbl.area_labels{vigrid};  end
        elseif isfield(opt,'viselec') && isfield(opt.viselec,'wang2015_atlas')
           vigrid = find(strcmp(opt.viselec.wang2015_atlas.elec_labels,gridList{ee}),1);
           if ~isempty(vigrid),  VA = opt.viselec.wang2015_atlas.area_labels{vigrid};  end
        end
        if strcmpi(opt.plot.addWangToTitle,'yes')&&~isempty(VA)&&~strcmp(VA,'none')
            areaList{ee} = sprintf('%s w:%s',areaList{ee},VA);
        end
        
        VA = [];
        if isTableCol(channels,'bensonarea')
           if iscell(channels.bensonarea),  VA = channels.bensonarea{igrid};
           else,                            VA = channels.bensonarea(igrid);
           end
        elseif isfield(opt,'viselec') && isfield(opt.viselec,'benson14_varea')
           vigrid = find(strcmp(opt.viselec.benson14_varea.elec_labels,gridList{ee}),1);
           if ~isempty(vigrid),  VA = opt.viselec.benson14_varea.area_labels{vigrid};  end
        end
        if strcmpi(opt.plot.addBensonToTitle,'yes')&&~isempty(VA)&&~strcmp(VA,'none')
            areaList{ee} = sprintf('%s b:%s',areaList{ee},VA);
            switch opt.plot.addEccToTitle
                case 'yes'
                    if isTableCol(channels,'bensoneccen')
                        VC = channels.bensoneccen(igrid);
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
if ~hasHDgrid
    v = reshape(1:FullGRID,[nCol,nRow])';
elseif ~opt.plot.RotGrid
    v = reshape(1:FullGRID,[nRow,nCol]);
else
    v = flipud(reshape(1:FullGRID,[nCol,nRow])'); 
end
inx = [];
for ee = 1:nFig
    inx{ee} = reshape(v(plRow*(ee-1)+1:min(nRow,plRow*(ee)),:)',1,[]);
end

%-- get gaussian properties of pRF & group updates
PRFs = cell(size(dat));
for ii = 1:ndats
    PRFs{ii} = [dat{ii}.ecc(:,1) .* cos(dat{ii}.ang(:,1)/180*pi), ...
                dat{ii}.ecc(:,1) .* sin(dat{ii}.ang(:,1)/180*pi), ...
                dat{ii}.rfsize(:,1)];
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
if ngroups>1
    lgndr = lgnds([1:length(lgnds)]>ngroups);
    lgnds = lgnds([1:length(lgnds)]<=ngroups);
else
    lgndr = lgnds(find(repidx)<=length(lgnds));
    lgnds = lgnds(find(~repidx)<=length(lgnds));
end

%-- set ylim
if ~isempty(opt.plot.XLim) && ~isempty(opt.plot.YLim)
        xmin = opt.plot.XLim(1) ./opt.plot.pix2deg;
        xmax = opt.plot.XLim(2) ./opt.plot.pix2deg;
        ymin = opt.plot.YLim(1) ./opt.plot.pix2deg;
        ymax = opt.plot.YLim(2) ./opt.plot.pix2deg;
else
    grididx = chanIdx(~isnan(chanIdx));
    xmin = prctile(prctile(PRFa(grididx,1,:)-max(opt.plot.sigma).*PRFa(grididx,3,:),5,3),5,1);
    xmax = prctile(prctile(PRFa(grididx,1,:)+max(opt.plot.sigma).*PRFa(grididx,3,:),95,3),95,1);
    ymin = prctile(prctile(PRFa(grididx,2,:)-max(opt.plot.sigma).*PRFa(grididx,3,:),5,3),5,1);
    ymax = prctile(prctile(PRFa(grididx,2,:)+max(opt.plot.sigma).*PRFa(grididx,3,:),95,3),95,1);
    if ~isempty(opt.plot.YLim)
        xmin = opt.plot.XLim(1) ./opt.plot.pix2deg;
        xmax = opt.plot.XLim(2) ./opt.plot.pix2deg;
        axlen = xmax - xmin;
        axcnt = mean([ymin,ymax]);
        ymin = round(axcnt,-round(log10(axlen)-2)) - axlen/2;
        ymax = round(axcnt,-round(log10(axlen)-2)) + axlen/2;
    elseif ~isempty(opt.plot.XLim)
        ymin = opt.plot.YLim(1) ./opt.plot.pix2deg;
        ymax = opt.plot.YLim(2) ./opt.plot.pix2deg;
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

%%% main
axmaxmin = max(xmax,ymax)-min(xmin,ymin);
axlen    = max(xmax-xmin,ymax-ymin);
axspd    = axmaxmin./axlen;
axplt    = linspace(min(xmin,ymin),max(xmax,ymax),ceil(opt.plot.resolution * axspd));
[xi, yi]=meshgrid(axplt,axplt);
for ee = 1:length(inx)
    %-- legend
    flglgnd = ~isempty(lgnds);
    %-- Decide how many subplots are needed
    nRow = opt.plot.nSubPlots(1);
    nCol = opt.plot.nSubPlots(2);
       
    %-- subplot
    hF = figure('Name', opt.plot.FigName); 
    set(hF, 'Position', [150 100 125*nCol 160*nRow]);%[150 100 2000 1250]
    for ichidx=1:length(inx{ee})
%       plotElectrodes = gridList(inx{ee}(ichidx));
      nameElectrodes = areaList(inx{ee}(ichidx));
      el = chanIdx(inx{ee}(ichidx));
      if ~isnan(el)
          %-- prepare gaussian image
          PRFexp = (bsxfun(@minus,xi,PRFs(el,1,:)).^2 + bsxfun(@minus,yi,PRFs(el,2,:)).^2)./(2.*PRFs(el,3,:).^2);
          PRFamp = 1;
          imPRFs = PRFamp  .* exp(-PRFexp);
          if opt.plot.isavg
              imPRFa = [];
              for ii=1:ngroups,    imPRFa = cat(3,imPRFa,nanmean(imPRFs(:,:,groups==ii),3));    end
              imPRFs = imPRFa;
          end
          PRFexp = (bsxfun(@minus,xi,PRFr(el,1,:)).^2 + bsxfun(@minus,yi,PRFr(el,2,:)).^2)./(2.*PRFr(el,3,:).^2);
          PRFamp = 1;
          imPRFr = PRFamp  .* exp(-PRFexp);
                    
          %-- Prepare axis       hA = subplot(nRow,nCol,ichidx);
          set(0,'CurrentFigure',hF);
          ax_sp = 1./([nCol nRow]).*[0.3 0.15];
          ax_siz = (1-ax_sp.*[1 2.0]) ./[nCol nRow];
          ax_pos = [ax_sp(1) + ax_siz(1)*(mod(ichidx-1,nCol)),...
                    ax_sp(2) + ax_siz(2)*(nRow-ceil(ichidx./nCol)),...
                    ax_siz(1)*0.70, ax_siz(2)*0.60];
          hA = subplot('Position',ax_pos);    daspect([1 1 1]);
          hold on
          
          %-- Plot
          hplts = [];
          switch opt.plot.Type
              case {'contour'}
                  plot([min(axplt) max(axplt)],[0 0],'k:');
                  plot([0 0],[min(axplt) max(axplt)],'k:');
                  drawellipse(0,0,0,opt.plot.resolution/2.*opt.plot.pix2deg,opt.plot.resolution/2.*opt.plot.pix2deg,[],[],'k:');
              case {'imagesc'}
                  ii=1;
                  imagesc(axplt.*opt.plot.pix2deg,axplt.*opt.plot.pix2deg,zeros(length(axplt)),...
                              opt.plot.options{groups(ii)}{:},[0 exp(-0.5.^2./2)]);
          end
          for ii = 1:size(imPRFs,3)
              switch opt.plot.Type
                  case {'contour'}
                      [~,hplts(ii)] = contour(axplt.*opt.plot.pix2deg,axplt.*opt.plot.pix2deg,imPRFs(:,:,ii),...
                                  exp(-opt.plot.sigma.^2./2), 'LineColor', plotCol(ii,:),'LineWidth',1.2,...
                                  opt.plot.options{groups(ii)}{:});
                      if ~opt.plot.isavg && (nanmax(nanmax(imPRFs(:,:,ii))) <= max(exp(-opt.plot.sigma.^2./2)))
                          plot(PRFs(el,1,ii).*opt.plot.pix2deg,PRFs(el,2,ii).*opt.plot.pix2deg,'x','LineWidth',1.5,'MarkerEdgeColor',plotCol(ii,:),'MarkerSize',8);
                      end
                  case {'imagesc'}
                      hplts(ii) = imagesc(axplt.*opt.plot.pix2deg,axplt.*opt.plot.pix2deg,imPRFs(:,:,ii),...
                                  'AlphaData',ones(size(imPRFs,[1 2]))./(size(imPRFs,3)+2*size(imPRFr,3)),...
                                  opt.plot.options{groups(ii)}{:},[0 exp(-0.5.^2./2)]);
              end
          end
          hpltr = [];
          for ii = 1:size(imPRFr,3)
              switch opt.plot.Type
                  case {'contour'}
                      [~,hpltr(ii)] = contour(axplt.*opt.plot.pix2deg,axplt.*opt.plot.pix2deg,imPRFr(:,:,ii),...
                                  exp(-opt.plot.sigma.^2./2),'LineColor', reptCol(groupr(ii),:),'LineWidth',2.4,...
                                  opt.plot.options{groupr(ii)}{:});
                      if ~opt.plot.isavg && (nanmax(nanmax(imPRFs(:,:,ii))) <= max(exp(-opt.plot.sigma.^2./2)))
                          plot(PRFr(el,1,ii).*opt.plot.pix2deg,PRFr(el,2,ii).*opt.plot.pix2deg,'x','LineWidth',2,'MarkerEdgeColor',reptCol(groupr(ii),:),'MarkerSize',15);
                      end
                  case {'imagesc'}
                      hpltr(ii) = imagesc(axplt.*opt.plot.pix2deg,axplt.*opt.plot.pix2deg,imPRFr(:,:,ii),...
                                  'AlphaData',ones(size(imPRFs,[1 2]))./(size(imPRFs,3)+2*size(imPRFr,3))*2,...
                                  opt.plot.rep_options{groupr(ii)}{:},[0 exp(-0.5.^2./2)]);
              end
          end
          axis xy equal
          
          switch opt.plot.Type
              case {'imagesc'}
                  plot([min(axplt) max(axplt)],[0 0],'w:');
                  plot([0 0],[min(axplt) max(axplt)],'w:');
                  drawellipse(0,0,0,opt.plot.resolution/2.*opt.plot.pix2deg,opt.plot.resolution/2.*opt.plot.pix2deg,[],[],'w:');
          end
          
          hold off
          
          %-- Set axis, title & legend
          set(hA,'XLim',[xmin xmax].*opt.plot.pix2deg,'YLim',[ymin ymax].*opt.plot.pix2deg,...
              'FontSize', opt.plot.fontSize);
          title(nameElectrodes);
%           setsubplotaxes(hA);
          if flglgnd
              lgndsidx = ~cellfun(@isempty,lgnds);      lgndridx = ~cellfun(@isempty,lgndr);
              hlg = legend([hplts(lgndsidx) hpltr(lgndridx)],[lgnds(lgndsidx) lgndr(lgndridx)],...
                  'FontSize',opt.plot.fontSize);
              hlg.NumColumns = nCol;
              hlg.Position([1,2]) = [ax_sp(1), hA.OuterPosition([2,4])*[1;1.02]];
              flglgnd = false;
          end
          drawnow limitrate nocallbacks;
      end
    end
    figlist{ee} = hF;
end
