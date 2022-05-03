% plotGridBarCommon
%   separate common script from ecog_plotGridBar functions

% 20220222 Yuasa

% Dependency: cellstrfind

%%
%-- distinguish HDgrid
hdGthresh = 64;

%-- correct elecnames (set '%02d' or '%03d')
numpltchan  = sum(cellfun(@(C) ~isempty(C), regexp(opt.channels.name,['^' whichGrid '\d+$'])));
[eleccat, elecnum] = strtok(opt.channels.name,int2str(0:9));
chanformat  = sprintf('%%s%%%02dd',max(cellfun(@length, elecnum)));
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
    gridList{ee} = sprintf(chanformat',whichGrid,ee);
    
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
              if eval(tickColChk)
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