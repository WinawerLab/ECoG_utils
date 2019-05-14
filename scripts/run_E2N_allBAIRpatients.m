

% close all;
if isempty(which('electrode_to_nearest_node'))
    tbUse('ECoG_utils');
end

%% run script for a specific patient (with plotmesh on) to look at electrodes on brain

%pIDs  = [648 645 661 668 674 682 692 708 718 723]; % BAIR patients
pIDs = 648;

% run script for each patient
specs = [];
out_all = cell(length(pIDs),1);
for i = 1:length(pIDs)
    
    specs.pID           = num2str(pIDs(i));
    specs.atlasNames    = [];%{'wang2015_atlas','wang15_mplbl','benson14_varea', 'benson14_eccen'}; 
                            % defaultclose al is all of these maps: 
                            % {'wang2015_atlas', ...
                            % 'benson14_varea', ...
                            % 'benson14_eccen', ...
                            % 'benson14_angle', ...
                            % 'benson14_sigma', ...
                            % 'template_areas'};
                            % NOTE: including benson14_varea is required to
                            % be able to obtain benson14_eccen, angle and sigma
    specs.plotmesh      = 'both';
    specs.plotelecs     = 'yes';
    specs.plotlabel     = 'no';  
    specs.patientPool   = 'BAIR';
    specs.elecFile      = sprintf('/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/preprocessed/sub-som%d/som%d_fsaverage_electrode_coordinates.tsv', pIDs(i), pIDs(i));
    specs.thresh = [];
    
    out_all{i} = electrode_to_nearest_node_fsaverage(specs);
end

