% plotGridSCCommon
%   separate common script from ecog_plotGridSC functions

% 20220222 Yuasa

% Dependency: cellstrfind

%%
%-- distinguish HDgrid
hdGthresh = 64;

%-- correct elecnames (set '%02d' or '%03d')
numpltchan  = sum(cellfun(@(C) ~isempty(C), regexp(opt.channels.name,['^' whichGrid '\d+$'])));
[eleccat, elecnum] = strtok(opt.channels.name,int2str(0:9));
if numpltchan > hdGthresh    % HDgrid
    chanformat = '%s%03d';
else
    chanformat = '%s%02d';
end
for el = 1:length(opt.channels.name)
    opt.channels.name{el} = sprintf(chanformat,eleccat{el},str2double(elecnum{el}));
end

%-- set grid parameters
if numpltchan > hdGthresh    % HDgrid
    nCol = 8; nRow = 16;
else
    nCol = 8; nRow = 8;
end
if opt.plot.RotGrid, tmp = nCol; nCol = nRow; nRow = tmp; end
FullGRID = nCol * nRow;

if isempty(opt.plot.nSubPlots)
    nFig = ceil(nRow / nCol ./ 2.2);
    opt.plot.nSubPlots = [ceil(nRow./nFig), nCol];
end
plRow  = opt.plot.nSubPlots(1);
nFig = ceil(nRow./plRow);

%-- Create a list of electrode names for the overall grid layout
gridList = cell(1,FullGRID);
grididx  = nan(1,FullGRID);
for ee = 1:FullGRID
    gridList{ee} = sprintf(chanformat,whichGrid,ee);
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
if ischar(opt.plot.CLim)
  switch opt.plot.CLim
    case 'maxmin',  yrange = [ymin ymax];
    case 'maxabs',  yrange = [-1 1].*max(abs([ymin ymax]));
    case 'minabs',  yrange = [-1 1].*min(abs([ymin ymax]));
    otherwise,      error('unknown parameter');
  end
else
    yrange = opt.plot.CLim;
end
setytick = opt.plot.colorBar && ...
    ~(ischar(opt.plot.CTick)&&strcmp(opt.plot.CTick,'auto'));
setyticklabel = opt.plot.colorBar && ...
    ~(ischar(opt.plot.CTickLabel{1})&&strcmp(opt.plot.CTickLabel{1},'auto'));
