function [channels] = bair_addVisualAtlasNamesToChannelTable(channels,visualelectrodes)

hasVisualElecs = 0;

% Create empty arrays to be added as columns to the channel table
wang         = repmat({'none'},[height(channels) 1]); 
benson       = repmat({'none'},[height(channels) 1]); 
hcp          = repmat({'none'},[height(channels) 1]); 
benson_eccen = nan([height(channels) 1]); 
benson_angle = nan([height(channels) 1]); 
benson_sigma = nan([height(channels) 1]);  
matched_node = nan([height(channels) 1]); 

% Check if there are any matched visual electrodes:
for atlas = {'wang2015_atlas','benson14_varea', 'wang15_mplbl','wang15_fplbl'}
	if ~isempty(visualelectrodes)
        if hasvalue(visualelectrodes,atlas{:})
            hasVisualElecs = 1;
        end
	end
end

if hasVisualElecs
    
    % Find all matched electrodes
    %-- Wang atlas
    viselec_wang = []; visarea_wang = []; matchednodes_wang = [];
    for wangname = {'wang2015_atlas', 'wang15_mplbl'}
        if hasvalue(visualelectrodes,wangname{:})
            viselec_wang = [viselec_wang visualelectrodes.(wangname{:}).elec_labels];
            visarea_wang = [visarea_wang visualelectrodes.(wangname{:}).area_labels];
            matchednodes_wang = [matchednodes_wang visualelectrodes.(wangname{:}).node_indices];
        end
    end

    %-- Wang probability    
    viselec_wang_prob = {}; visarea_wang_prob = {}; visval_wang_prob = {}; matchednodes_wang_prob = {};
    if hasvalue(visualelectrodes,'wang15_fplbl')
        visarea_wang_prob = [visarea_wang_prob visualelectrodes.wang15_fplbl.area_names];
        viselec_wang_prob = [viselec_wang_prob visualelectrodes.wang15_fplbl.elec_labels];
        visval_wang_prob  = [visval_wang_prob visualelectrodes.wang15_fplbl.node_values];
        matchednodes_wang_prob = [matchednodes_wang_prob visualelectrodes.wang15_fplbl.node_indices];
    end
    wang_prob = zeros([height(channels) length(visarea_wang_prob)]);
    
    %-- Benson atlas
    viselec_benson = []; 
    if hasvalue(visualelectrodes,'benson14_varea')
        viselec_benson = [viselec_benson visualelectrodes.benson14_varea.elec_labels];
    end
    
    %-- Benson HCP atlas
    viselec_hcp = []; 
    if hasvalue(visualelectrodes,'benson20_mplbl')
        viselec_hcp = [viselec_hcp visualelectrodes.benson20_mplbl.elec_labels];
    end
    
    %-- Benson HCP probability
    viselec_hcp_prob = {}; visarea_hcp_prob = {}; visval_hcp_prob = {}; matchednodes_hcp_prob = {};
    if hasvalue(visualelectrodes,'benson20_fplbl')
        visarea_hcp_prob = [visarea_hcp_prob visualelectrodes.benson20_fplbl.area_names];
        viselec_hcp_prob = [viselec_hcp_prob visualelectrodes.benson20_fplbl.elec_labels];
        visval_hcp_prob  = [visval_hcp_prob visualelectrodes.benson20_fplbl.node_values];
        matchednodes_hcp_prob = [matchednodes_hcp_prob visualelectrodes.benson20_fplbl.node_indices];
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

fprintf('[%s] Adding visual matches to channel table\n',mfilename);

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

end

function tf = hasvalue(S,field)
tf = isfield(S,field) && ~isempty(S.(field));
end
