function figlist = ecog_plotGridPRF(whichGrid, opt, varargin)

% p = ECOG_PLOTGRIDPRF(whichGrid, option, result1)
% p = ECOG_PLOTGRIDPRF(whichGrid, option, result1, result2,...)
% p = ECOG_PLOTGRIDPRF(whichElectrodes, option, result1)
% ECOG_PLOTGRIDPRF draw pRF in the entire grid.
% 
% <whichGrid> can be 'G', 'GA' or 'GB'.
% Otherwise it wll be interpreted as electrodes names (some options will be ignored).
% 
% <result> is a structure of pRF information, requiring: 
%   result.ecc:    Nx1 array with pRF eccentricity
%   result.ang:    Nx1 array with pRF angle [degree: 0deg is on positive x-axis and increase anti-clockwise]
%   result.rfsize: Nx1 array with pRF size
% each result can be Mx1 cell-array with pRF structures, or Nx1xM array
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
%       addCenter        = 'yes' or 'no'(default)  % show center as '.' if plot type is 'contour'
% 
%       addChsToTitle    = 'yes'(default) or 'no'
%       addWangToTitle   = 'yes'(default) or 'no'
%       addBensonToTitle = 'yes'(default) or 'no'
%       addR2ToTitle     = 'yes' or 'no'(default)  % show R2 or xval of result1
%       addEccToTitle    = 'yes' or 'no'(default)  % show eccentricity of result1
%       addSbjToTitle    = 'yes' or 'no'(default)  % show subject name
%       nSubPlots        = [m n], plot m-by-n grid electrodes in a figure.
%       RotGrid          = true or false (default), if true, electrodes are
%                          aligned along row instead of column.
%       labelCol         = color or cell-array of colors, color of the titles
%                          in each axis. defualt is black.
%       FigName          = string, name of the figure.
%       FigSizRate       = scalar (default=1).
%       fontSize         = scalar, font size in the figure.
%       legend           = string or cell-array of strings, legend labels.
%       options          = any other options for each axis.
%   option.plotbenson    = 'yes' or 'no'(default)  % plot results estimated from benson atlas
%   option.subjet        = string, subject name is helpful to interpret benson prf.
% 
% Followings are unexplained options for option.plot: 
%       representative
%       rep_options
%       showaxis

% option should include 'channels'

% Dependency: <analyzePRF>, plotGridCommon, SetDefault, cellstrfind, arrangeinrect, istablefield
% Develop: MATLAB2020b

% 20191106 Yuasa
% 20200206 Yuasa: change arguments structures (not compatible to prior version)
% 20200309 Yuasa: enable to plot no HD grid data
% 20200331 Yuasa: add 'plotbenson' option
% 20200925 Yuasa: solve problem using strings in channels
% 20201110 Yuasa: bug fix
% 20201228 Yuasa: update for data having multiple pRF results
% 20210316 Yuasa: minor change
% 20210419 Yuasa: use tiledlayout
% 20210701 Yuasa: use plotGridCommon
% 20220307 Yuasa: remove patient dependent code

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
SetDefault('opt.plot.addCenter', 'no', 0);
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
SetDefault('opt.plot.FigSizRate', 1, 0);
SetDefault('opt.plot.RotGrid', true, 0);
SetDefault('opt.plot.colors', [], 1);
SetDefault('opt.plot.labelCol', [], 1);
SetDefault('opt.plot.representative', [], 1);
SetDefault('opt.plot.legend', {}, 1,'cell');
SetDefault('opt.plot.options', {}, 1);
SetDefault('opt.plot.rep_options', {}, 1);
SetDefault('opt.plot.showaxis', true, 0);

SetDefault('whichGrid','*',0);

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

%-- data correction (align for multi group)
ngroups = nargin-2;
groups  = [];
for ii=1:ngroups
    if iscell(varargin{ii}),  varargin{ii} = reshape(varargin{ii},[],1);
    else,                     varargin{ii} = boot2cell(varargin{ii});
    end
    groups  = [groups; ii*ones(size(varargin{ii}))];
end
dat      = cat(1,varargin{:});
% for ii=1:numel(dat)
    
ndats    = numel(dat);

%-- Benson reference
if plotbenson
    assert(istablefield(channels,'hemisphere'),'hemisphere column must be included in channel table');
    islefths = ismember(channels.hemisphere(elec),'L');     % right visual field
        
    ngroups = ngroups+1;
    groups  = [groups; ngroups];
    ndats   = ndats+1;
    
    %-- Correct degree to pixel, Benson angle to analyzePRF angle
    dat{ndats}.ecc      = channels.bensoneccen./opt.plot.pix2deg;
    dat{ndats}.ang      = mod(((-1).^islefths).*channels.bensonangle+90,360);
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

%%% run plotGridCommon
plotGridCommon;

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
       
    %-- figure  % row+30 for legend  % -20*row if no-axis
    hF = figure('Name', opt.plot.FigName); 
    hF.Position = [150 100 125*nCol (160-20.*~opt.plot.showaxis)*nRow+30] ...
                    .* [1 1 opt.plot.FigSizRate(1).^1.38 opt.plot.FigSizRate(1)];
    hT = tiledlayout(nRow,nCol,'TileSpacing','compact','Padding','compact');
    %-- subplot
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
              imPRFa = []; parPRFs = [];
              for ii=1:ngroups
                  imPRFa  = cat(3,imPRFa,nanmean(imPRFs(:,:,groups==ii),3));
                  parPRFs = cat(3,parPRFs,nanmean(PRFs(el,:,groups==ii),3));
              end
              imPRFs = imPRFa;
          else
              parPRFs = PRFs(el,:,:);
          end
          PRFexp = (bsxfun(@minus,xi,PRFr(el,1,:)).^2 + bsxfun(@minus,yi,PRFr(el,2,:)).^2)./(2.*PRFr(el,3,:).^2);
          PRFamp = 1;
          imPRFr = PRFamp  .* exp(-PRFexp);
          parPRFr = PRFr(el,:,:);
                    
          %-- Prepare axis       hA = subplot(nRow,nCol,ichidx);
          set(0,'CurrentFigure',hF);
          ax_sp = 1./([nCol nRow]).*[0.3 0.15];
%           ax_siz = (1-ax_sp.*[1 2.0]) ./[nCol nRow];
%           ax_pos = [ax_sp(1) + ax_siz(1)*(mod(ichidx-1,nCol)),...
%                     ax_sp(2) + ax_siz(2)*(nRow-ceil(ichidx./nCol)),...
%                     ax_siz(1)*0.70, ax_siz(2)*0.60];
%           hA = subplot('Position',ax_pos);    daspect([1 1 1]);
% %           hA = subplot_er(nRow,nCol,ichidx);    daspect([1 1 1]);
          hA = nexttile(ichidx);    daspect([1 1 1]);
          hold on
          
          %-- Plot
          %%-- plot all pRFs except for representative pRFs
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
                      if (nanmax(nanmax(imPRFs(:,:,ii))) <= max(exp(-opt.plot.sigma.^2./2)))
                          plot(parPRFs(1,1,ii).*opt.plot.pix2deg,parPRFs(1,2,ii).*opt.plot.pix2deg,'x','LineWidth',1.5,'MarkerEdgeColor',plotCol(ii,:),'MarkerSize',8);
                      elseif strcmpi(opt.plot.addCenter,'yes')
                          plot(parPRFs(1,1,ii).*opt.plot.pix2deg,parPRFs(1,2,ii).*opt.plot.pix2deg,'.','LineWidth',1.5,'MarkerEdgeColor',plotCol(ii,:),'MarkerSize',8);
                      end
                  case {'imagesc'}
                      hplts(ii) = imagesc(axplt.*opt.plot.pix2deg,axplt.*opt.plot.pix2deg,imPRFs(:,:,ii),...
                                  'AlphaData',ones(size(imPRFs,[1 2]))./(size(imPRFs,3)+2*size(imPRFr,3)),...
                                  opt.plot.options{groups(ii)}{:},[0 exp(-0.5.^2./2)]);
              end
          end
          %%-- plot representative pRFs
          hpltr = [];
          for ii = 1:size(imPRFr,3)
              switch opt.plot.Type
                  case {'contour'}
                      [~,hpltr(ii)] = contour(axplt.*opt.plot.pix2deg,axplt.*opt.plot.pix2deg,imPRFr(:,:,ii),...
                                  exp(-opt.plot.sigma.^2./2),'LineColor', reptCol(groupr(ii),:),'LineWidth',2.4,...
                                  opt.plot.rep_options{groupr(ii)}{:});
                      if (nanmax(nanmax(imPRFs(:,:,ii))) <= max(exp(-opt.plot.sigma.^2./2)))
                          plot(parPRFr(1,1,ii).*opt.plot.pix2deg,parPRFr(1,2,ii).*opt.plot.pix2deg,'x','LineWidth',2,'MarkerEdgeColor',reptCol(groupr(ii),:),'MarkerSize',16);
                      elseif strcmpi(opt.plot.addCenter,'yes')
                          plot(parPRFr(1,1,ii).*opt.plot.pix2deg,parPRFr(1,2,ii).*opt.plot.pix2deg,'.','LineWidth',2,'MarkerEdgeColor',reptCol(groupr(ii),:),'MarkerSize',16);
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
%               hlg.Position([1,2]) = [ax_sp(1), hA.OuterPosition([2,4])*[1;1.02]];
%               hlg.Location='northoutside';
              hlg.Layout.Tile = 'North';
              flglgnd = false;
          end
          if ~opt.plot.showaxis
              hA.XTick = [];   hA.YTick = [];   hA.ZTick = [];
              hA.XAxis.Visible = 'off';
              hA.YAxis.Visible = 'off';
              hA.ZAxis.Visible = 'off';
          end
          drawnow limitrate nocallbacks;
      end
    end
    figlist{ee} = hF;
end

function results_cell = boot2cell(results)
prfflds = {'ang','ecc','expt','rfsize','R2','gain','testperformance','aggregatedtestperformance','xR2','xval','resnorms','numiters','meanvol','params'};
lstfids = reshape(fieldnames(results),1,[]);
chkflds = lstfids(ismember(lstfids,prfflds));
for fld = chkflds
    if isempty(results.(fld{:})),  chkflds(ismember(chkflds,fld{:})) = []; end
end

%-- check iteration
if isempty(chkflds), 	niter = 1;
else,               	niter = size(results.(chkflds{1}),ndims(results.(chkflds{1})));
end

results_cell = cell(niter,1);
for ii = 1:niter
    %-- copy other fields
    results_cell{ii} = rmfield(results,chkflds);
    %-- copy prf fields
    for fld = chkflds
        fldidx = sprintf('%s%d',repmat(':,',1,ndims(results.(fld{:}))-1),ii);
        results_cell{ii}.(fld{:}) = eval(sprintf('%s(%s)','results.(fld{:})',fldidx));
    end
end

