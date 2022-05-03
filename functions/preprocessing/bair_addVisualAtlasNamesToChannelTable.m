function [channels,chanidx] = bair_addVisualAtlasNamesToChannelTable(channels,visualelectrodes)
% [channels,index] = bair_addVisualAtlasNamesToChannelTable(channels,electrodes)
%   add atlas information to channels
% 
% Input
%     channels:  	channel table
%     electrodes: 	electrode table including atlas information
%                       or
%                   data structure including electrode matching information 
%                   (compatible modefor outputs of electrode_to_nearest_node)
% 
% Output
%     channels:  	channel table with atlas
%                   if electrodes is table, then only channels with
%                   coordinates are output
%     index:        index to match input and output channels
% 
% See also, bidsEcogMatchElectrodesToAtlas
% 
% IG, BAIR 2019; K.Yuasa, 2020-22

% Atlas names
wangatlas   = {'wang2015_atlas', '*wang15_mplbl','*wang15_fplbl_norm*'};
wangprobs   = {'*wang15_fplbl'};
bensonatlas = {'*benson14_varea'};
hcpatlas    = {'benson20_mplbl','benson20_fplbl_norm*'};
hcpprobs    = {'benson20_fplbl'};
bensonprfs  = {'benson14_eccen','benson14_angle','benson14_sigma'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  if visualelectrodes is table   %%%
if istable(visualelectrodes)
    fprintf('[%s] Adding atlas information to channel table\n',mfilename);
    
    % Solve fplbl atlas names
    %-- Wang probability
    wangprobs = [wangprobs, cellfun(@(C) sprintf('%s_*',C),wangprobs,'UniformOutput',false)];
    wangprobs = [wangprobs(1), setdiff(chkfields(visualelectrodes,wangprobs),...
                                       chkfields(visualelectrodes,wangatlas),'stable')];
    %-- HCP probability
    hcpprobs  = [hcpprobs, cellfun(@(C) sprintf('%s_*',C),hcpprobs,'UniformOutput',false)];
    hcpprobs  = [hcpprobs(1), setdiff(chkfields(visualelectrodes,hcpprobs),...
                                      chkfields(visualelectrodes,hcpatlas),'stable')];
    
    % Select matched channels
    [~, elecidx, chanidx] = intersect(visualelectrodes.name, channels.name,'stable');
    channels    = channels(chanidx,:);
    visualelectrodes  = visualelectrodes(elecidx,:);
    
    % Get valid atlas names
    atlas = [wangatlas,wangprobs,bensonatlas,hcpatlas,hcpprobs,bensonprfs];
    atlasfields = chkfields(visualelectrodes,atlas);
    
    % Get additional fields from electrods
    chanfields = tbfieldnames(channels);
    elecfields = tbfieldnames(visualelectrodes);
    [extrfields] = setdiff(elecfields,[chanfields,atlasfields],'stable');
    
    % Add addtional info to channels
    channels = horzcat(channels, visualelectrodes(:,extrfields));
    
    % Add electrode atlas info to channels
    %-- Benson atlas
    for atlasname = chkfields(visualelectrodes,bensonatlas)
        if hasvalue(visualelectrodes,atlasname{:})
            channels.bensonarea  = visualelectrodes.(atlasname{:});
        end
    end
    
    %-- Wang atlas
    for atlasname = chkfields(visualelectrodes,wangatlas)
        if hasvalue(visualelectrodes,atlasname{:})
            channels.wangarea  = visualelectrodes.(atlasname{:});
        end
    end
    
    %-- Benson HCP atlas
    for atlasname = chkfields(visualelectrodes,hcpatlas)
        if hasvalue(visualelectrodes,atlasname{:})
            channels.hcparea  = visualelectrodes.(atlasname{:});
        end
    end
    
    %-- Wang probability
    for atlasname = chkfields(visualelectrodes,wangprobs)
        if hasvalue(visualelectrodes,atlasname{:})
            visarea  = regexprep(atlasname{:},regexptranslate('wildcard',[wangprobs{1} '_']),'');
            channels.(['wangprob_' visarea])  = visualelectrodes.(atlasname{:});
        end
    end
    
    %-- Benson HCP probability
    for atlasname = chkfields(visualelectrodes,hcpprobs)
        if hasvalue(visualelectrodes,atlasname{:})
            visarea  = regexprep(atlasname{:},regexptranslate('wildcard',[hcpprobs{1} '_']),'');
            channels.(['hcpprob_' visarea])  = visualelectrodes.(atlasname{:});
        end
    end
    
    %-- Benson pRFs
    for atlasname = chkfields(visualelectrodes,bensonprfs)
        if hasvalue(visualelectrodes,atlasname{:})
            if endsWith(atlasname{:},'_eccen')
                channels.bensoneccen  = visualelectrodes.(atlasname{:});
            elseif endsWith(atlasname{:},'_angle')
                channels.bensonangle  = visualelectrodes.(atlasname{:});
            elseif endsWith(atlasname{:},'_sigma')
                channels.bensonsigma  = visualelectrodes.(atlasname{:});
            end
        end
    end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  if visualelectrodes is struct  %%%
else
        
    % Create empty arrays to be added as columns to the channel table
    wang         = repmat({'none'},[height(channels) 1]);
    benson       = repmat({'none'},[height(channels) 1]);
    hcp          = repmat({'none'},[height(channels) 1]);
    benson_eccen = nan([height(channels) 1]);
    benson_angle = nan([height(channels) 1]);
    benson_sigma = nan([height(channels) 1]);
    matched_node = nan([height(channels) 1]);
    
    hasVisualElecs = 0;
    % Check if there are any matched visual electrodes:
    atlas = [wangatlas,wangprobs,bensonatlas,hcpatlas,hcpprobs];
    for atlasname = chkfields(visualelectrodes,atlas)
        if ~isempty(visualelectrodes)
            if hasvalue(visualelectrodes,atlasname{:})
                hasVisualElecs = 1;
            end
        end
    end

    if hasVisualElecs
        
        fprintf('[%s] Adding visual matches to channel table\n',mfilename);

        % Find all matched electrodes
        %-- Wang atlas
        viselec_wang = []; visarea_wang = []; matchednodes_wang = [];    
        for atlasname = chkfields(visualelectrodes,wangatlas)
            if hasvalue(visualelectrodes,atlasname{:})
                viselec_wang = [viselec_wang visualelectrodes.(atlasname{:}).elec_labels];
                visarea_wang = [visarea_wang visualelectrodes.(atlasname{:}).area_labels];
                matchednodes_wang = [matchednodes_wang visualelectrodes.(atlasname{:}).node_indices];
            end
        end

        %-- Wang probability    
        viselec_wang_prob = {}; visarea_wang_prob = {}; visval_wang_prob = {}; matchednodes_wang_prob = {};
        for atlasname = chkfields(visualelectrodes,wangprobs)
            if hasvalue(visualelectrodes,atlasname{:})
                visarea_wang_prob = [visarea_wang_prob visualelectrodes.wang15_fplbl.area_names];
                viselec_wang_prob = [viselec_wang_prob visualelectrodes.wang15_fplbl.elec_labels];
                visval_wang_prob  = [visval_wang_prob visualelectrodes.wang15_fplbl.node_values];
                matchednodes_wang_prob = [matchednodes_wang_prob visualelectrodes.wang15_fplbl.node_indices];
            end
        end
        wang_prob = zeros([height(channels) length(visarea_wang_prob)]);

        %-- Benson atlas
        viselec_benson = []; 
        for atlasname = chkfields(visualelectrodes,bensonatlas)
            if hasvalue(visualelectrodes,atlasname{:})
                viselec_benson = [viselec_benson visualelectrodes.benson14_varea.elec_labels];
            end
        end

        %-- Benson HCP atlas
        viselec_hcp = []; 
        for atlasname = chkfields(visualelectrodes,hcpatlas)
            if hasvalue(visualelectrodes,atlasname{:})
                viselec_hcp = [viselec_hcp visualelectrodes.benson20_mplbl.elec_labels];
            end
        end

        %-- Benson HCP probability
        viselec_hcp_prob = {}; visarea_hcp_prob = {}; visval_hcp_prob = {}; matchednodes_hcp_prob = {};
        for atlasname = chkfields(visualelectrodes,hcpprobs)
            if hasvalue(visualelectrodes,atlasname{:})
                visarea_hcp_prob = [visarea_hcp_prob visualelectrodes.benson20_fplbl.area_names];
                viselec_hcp_prob = [viselec_hcp_prob visualelectrodes.benson20_fplbl.elec_labels];
                visval_hcp_prob  = [visval_hcp_prob visualelectrodes.benson20_fplbl.node_values];
                matchednodes_hcp_prob = [matchednodes_hcp_prob visualelectrodes.benson20_fplbl.node_indices];
            end
        end
        hcp_prob = zeros([height(channels) length(visarea_hcp_prob)]);

        % Get visual elec indices
        [viselec] = unique([viselec_wang viselec_benson viselec_wang_prob{:} viselec_hcp visarea_hcp_prob{:}])';

        for ii = 1:length(viselec)

            elInx = ecog_matchChannels(viselec{ii}, channels.name)';
            if elInx>0
                %name{ii} = viselec{ii};
                %-- Wang atlas
                wangInx = find(strcmp(viselec_wang, viselec{ii}));
                if ~isempty(wangInx)
                    wang{elInx} = visarea_wang{wangInx(1)};
                    matched_node(elInx) = matchednodes_wang(wangInx(1));           
                end
                %-- Benson atlas
                bensonInx = find(strcmp(viselec_benson, viselec{ii}));
                if ~isempty(bensonInx)
                    benson{elInx} = visualelectrodes.benson14_varea.area_labels{bensonInx};
                    benson_eccen(elInx) = visualelectrodes.benson14_varea.node_eccen(bensonInx);
                    benson_angle(elInx) = visualelectrodes.benson14_varea.node_angle(bensonInx);
                    benson_sigma(elInx) = visualelectrodes.benson14_varea.node_sigma(bensonInx);
                    matched_node(elInx) = visualelectrodes.benson14_varea.node_indices(bensonInx);
                end
                %-- Wang probability
                for jj = 1:length(visarea_wang_prob)
                    wangInx = find(strcmp(viselec_wang_prob{jj}, viselec{ii}));
                    if ~isempty(wangInx)
                        wang_prob(elInx,jj) = visval_wang_prob{jj}(wangInx(1));
                        matched_node(elInx) = matchednodes_wang_prob{jj}(wangInx(1));
                    end
                end
                %-- Benson HCP atlas
                bensonInx = find(strcmp(viselec_hcp, viselec{ii}));
                if ~isempty(bensonInx)
                    hcp{elInx} = visualelectrodes.benson20_mplbl.area_labels{bensonInx};
                    matched_node(elInx) = visualelectrodes.benson20_mplbl.node_indices(bensonInx);
                end
                %-- Benson HCP probability
                for jj = 1:length(visarea_hcp_prob)
                    bensonInx = find(strcmp(viselec_hcp_prob{jj}, viselec{ii}));
                    if ~isempty(bensonInx)
                        hcp_prob(elInx,jj) = visval_hcp_prob{jj}(bensonInx(1));
                        matched_node(elInx) = matchednodes_hcp_prob{jj}(bensonInx(1));
                    end
                end
            end
        end
    end

    % Add area matches to channel table
    channels.bensonarea  = benson;
    channels.wangarea    = wang;
    channels.hcparea     = hcp;
    channels.bensoneccen = benson_eccen;
    channels.bensonangle = benson_angle;
    channels.bensonsigma = benson_sigma;
    channels.matchednode = matched_node;
    for jj = 1:length(visarea_wang_prob)
        channels.(['wangprob_' visarea_wang_prob{jj}]) = wang_prob(:,jj);
    end
    for jj = 1:length(visarea_hcp_prob)
        channels.(['hcpprob_' visarea_hcp_prob{jj}]) = hcp_prob(:,jj);
    end

    chanidx = [];

    
end
end

function [fields] = chkfields(S,fields)
if istable(S), S = summary(S); end
if ~iscell(fields), fields = {fields}; end
newfileds = {};
infld = 1;
for ifld = 1:length(fields)
    field = fields{ifld};
    fieldparts = strsplit(field,'*');
    if length(fieldparts) == 1
        if isfield(S,field)
            newfileds{infld} = field;
            infld = infld + 1;
        end
    else
        fnames = reshape(fieldnames(S),[],1);
        fidx   = true(size(fnames));
        for ifps = 1:length(fieldparts)
            if  ~isempty(fieldparts{ifps})
                if ifps == 1
                    fidx = [fidx, startsWith(fnames,fieldparts{ifps})];
                elseif ifps == length(fieldparts)
                    fidx = [fidx, endsWith(fnames,fieldparts{ifps})];
                else
                    fidx = [fidx, contains(fnames,fieldparts{ifps})];
                end
            end
        end
        fidx = reshape(find(all(fidx,2)),1,[]);
        for ifname = fidx
            newfileds{infld} = fnames{ifname};
            infld = infld + 1;
        end
    end
end
fields = reshape(newfileds,1,[]);
end

function tf = hasvalue(S,field)
if istable(S), S = summary(S); end
tf = isfield(S,field) && ~isempty(S.(field));
end

function [fields] = tbfieldnames(T)
fields = fieldnames(summary(T));
fields = reshape(fields,1,[]);
end
