% plotGridBarSimpleCommon
%   separate common script from ecog_plotGridSpectra, ecog_plotGridTimecourses

% 20220222 Yuasa

%%
%-- distinguish HDgrid
hdGthresh = 64;

%-- get electrodes list
chanIdx     = cellfun(@(C) ~isempty(C), regexp(channels.name,['^' whichGrid '\d+$']));
numpltchan  = sum(chanIdx);

%-- correct elecnames (set '%03d')
[eleccat, elecnum] = strtok(channels.name,int2str(0:9));
chanformat  = sprintf('%%s%%%02dd',max(cellfun(@length, elecnum)));
for el = 1:length(channels.name)
    channels.name{el} = sprintf(chanformat,eleccat{el},str2double(elecnum{el}));
end
whichElectrodes = channels.name(chanIdx);

%-- set grid parameters
if numpltchan > hdGthresh    % HDgrid
    nCol = 8; nRow = 16;
else
    nCol = 8; nRow = 8;
end
if specs.plot.RotGrid, tmp = nCol; nCol = nRow; nRow = tmp; end
FullGRID = nCol * nRow;
if isempty(specs.plot.nSubPlots)
    nFig = ceil(nRow / nCol ./ 1.2);
    specs.plot.nSubPlots = [ceil(nRow./nFig), nCol];
end
plRow  = specs.plot.nSubPlots(1);
nFig = ceil(nRow./plRow);

%-- Create a list of electrode names for the overall grid layout
gridList = [];
for ee = 1:FullGRID
    chanName = sprintf(chanformat,whichGrid,ee);
        
    if any(strcmp(chanName,whichElectrodes))
        gridList{ee} = chanName;
    else
        gridList{ee} = [chanName '-nodata'];
    end
end
 
%-- Plot response per electrode, in two separate figures (top and bottom)
v = reshape(1:FullGRID,[nRow,nCol]);
inx = [];
for ee = 1:nFig
    inx{ee} = reshape(v(plRow*(ee-1)+1:min(nRow,plRow*(ee)),:)',1,[]);
end
