% plotGridCommon
%   separate common script from ecog_plotGrid functions

% 20210630 Yuasa
% 20220222 Yuasa: reorganize grid and HDgrid

% SetDefault('opt.plot.addChsToTitle', 'yes', 0);
% SetDefault('opt.plot.addWangToTitle', 'yes', 0);
% SetDefault('opt.plot.addBensonToTitle', 'yes', 0);
% SetDefault('opt.plot.addR2ToTitle', 'no', 0);
% SetDefault('opt.plot.addEccToTitle', 'no', 0);
% SetDefault('opt.plot.addSbjToTitle', 'no', 0);

% Dependency: cellstrfind

%%
%-- distinguish HDgrid
hdGthresh = 64;

%-- check grid
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

%-- correct elecnames (set '%02d' or '%03d')
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
    if isempty(opt.plot.nSubPlots)
        [nCol,nRow] = arrangeinrect(numel(gridList),2,[4,1]);
    elseif opt.plot.nSubPlots(2) == 0
        nRow = opt.plot.nSubPlots(1);
        nCol = ceil(numel(gridList)./opt.plot.nSubPlots(1));
    elseif prod(opt.plot.nSubPlots) < numel(gridList)
        nRow = ceil(numel(gridList)./opt.plot.nSubPlots(2));
        nCol = opt.plot.nSubPlots(2);
    else
        nRow = opt.plot.nSubPlots(1);
        nCol = opt.plot.nSubPlots(2);
    end
end
if hasGrid && opt.plot.RotGrid, tmp = nCol; nCol = nRow; nRow = tmp; end
FullGRID = nCol * nRow;

if isempty(opt.plot.nSubPlots)
    nFig = ceil(nRow / nCol ./ 1.2);
    opt.plot.nSubPlots = [ceil(nRow./nFig), nCol];
elseif opt.plot.nSubPlots(1) == 0
    opt.plot.nSubPlots(1) = nRow;
elseif opt.plot.nSubPlots(2) == 0
    opt.plot.nSubPlots(2) = nCol;
end
plRow  = opt.plot.nSubPlots(1);
nFig = ceil(nRow./plRow);

%-- Check cross-validation
if exist('result','var')
    prfdat = result;
elseif exist('dat','var') && iscell(dat)
    prfdat = dat{1};
else
    prfdat = option;
end
if isfield(prfdat,'xval') && ~all(isnan(prfdat.xval))
    R2field = 'xval';
elseif isfield(prfdat,'aggregatedtestperformance') && ~all(isnan(prfdat.aggregatedtestperformance))
    R2field = 'aggregatedtestperformance';
elseif isfield(prfdat,'R2')
    R2field = 'R2';
else
    opt.plot.addR2ToTitle = 'no';
end

if strcmpi(opt.plot.addSbjToTitle,'yes')
    subjects = cellstr(channels.subject_name);
end

%-- Create a list of electrode names for the overall grid layout
areaList = [];
isTableCol = @(t, thisCol) ismember(thisCol, t.Properties.VariableNames);
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
    
    %-- get Visual Area Label
    igrid = chanIdx(ee);
    if strcmpi(opt.plot.addR2ToTitle,'yes') && ~isnan(igrid)
        if strcmpi(opt.plot.addChsToTitle,'yes')
            areaList{ee} = sprintf('%s (%.1f%%)',gridList{ee},prfdat.(R2field)(igrid));
        else
            areaList{ee} = sprintf('(%.1f%%)',prfdat.(R2field)(igrid));
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
