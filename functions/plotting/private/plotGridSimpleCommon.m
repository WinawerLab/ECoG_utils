% plotGridBarSimpleCommon
%   separate common script from ecog_plotGridSpectra, ecog_plotGridTimecourses

% 20220222 Yuasa

%%
%-- distinguish HDgrid
hdGthresh = 64;

%-- check grid & get electrodes list
if ischar(whichGrid) && ismember(whichGrid,{'G','GA','GB'})
    hasGrid   = true;
    gridList    = [];
    chanIdx     = [];
    numpltchan  = sum(cellfun(@(C) ~isempty(C), regexp(channels.name,['^' whichGrid '\d+$'])));
else
    hasGrid   = false;
    if islogical(whichGrid)
      chanIdx     = find(whichGrid);
    elseif isnumeric(whichGrid) && isequaln(whichGrid,round(whichGrid))
      chanIdx     = whichGrid;
    else    % list of channel names
      chanIdx     = cellstrfind(channels.name, whichGrid,0);
    end
    gridList    = channels.name(chanIdx);
    whichGrid = 'Channels';
    numpltchan  = numel(chanIdx);
end

%-- correct elecnames (set '%03d')
if hasGrid
[eleccat, elecnum] = strtok(channels.name,int2str(0:9));
chanformat  = sprintf('%%s%%%02dd',max(cellfun(@length, elecnum)));
for el = 1:length(channels.name)
    channels.name{el} = sprintf(chanformat,eleccat{el},str2double(elecnum{el}));
end
end

%-- set grid parameters
if hasGrid
    if numpltchan > hdGthresh    % HDgrid
        nCol = 8; nRow = 16;
    else
        nCol = 8; nRow = 8;
    end
else
    if isempty(specs.plot.nSubPlots)
        [nCol,nRow] = arrangeinrect(numel(gridList),2,[4,1]);
    elseif specs.plot.nSubPlots(2) == 0
        nRow = specs.plot.nSubPlots(1);
        nCol = ceil(numel(gridList)./specs.plot.nSubPlots(1));
    elseif prod(specs.plot.nSubPlots) < numel(gridList)
        nRow = ceil(numel(gridList)./specs.plot.nSubPlots(2));
        nCol = specs.plot.nSubPlots(2);
    else
        nRow = specs.plot.nSubPlots(1);
        nCol = specs.plot.nSubPlots(2);
    end
end
if hasGrid && specs.plot.RotGrid, tmp = nCol; nCol = nRow; nRow = tmp; end
FullGRID = nCol * nRow;

if isempty(specs.plot.nSubPlots)
    nFig = ceil(nRow / nCol ./ 1.2);
    specs.plot.nSubPlots = [ceil(nRow./nFig), nCol];
elseif specs.plot.nSubPlots(1) == 0
    specs.plot.nSubPlots(1) = nRow;
elseif specs.plot.nSubPlots(2) == 0
    specs.plot.nSubPlots(2) = nCol;
end
plRow  = specs.plot.nSubPlots(1);
nFig = ceil(nRow./plRow);

%-- Create a list of electrode names for the overall grid layout
% gridList = [];
for ee = 1:FullGRID
    if hasGrid
        
      gridList{ee} = sprintf(chanformat,whichGrid,ee);
      
      igrid = find(strcmp(channels.name,gridList{ee}),1);
      if isempty(igrid),    chanIdx(ee)  = nan;
      else,                 chanIdx(ee)  = igrid;
      end
      
    elseif ee > numel(gridList)
        
      gridList{ee} = 'none';
      chanIdx(ee)  = nan;
      
    end
end
 
%-- Plot response per electrode, in two separate figures (top and bottom)
if ~hasGrid
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
