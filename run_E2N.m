% close all;
tbUse('ECoG_utils');

%% run script for a specific patient (with plotmesh on) to look at electrodes on brain

specs = [];

specs.pID           = '668'; % patient ID number
%specs.atlasNames    = {'wang2015_atlas', 'benson14_varea', 'benson14_eccen'}; 
                        % default is all of these maps: 
                        % {'wang2015_atlas', ...
                        % 'benson14_varea', ...
                        % 'benson14_eccen', ...
                        % 'benson14_angle', ...
                        % 'benson14_sigma', ...
                        % 'template_areas'};
                        % NOTE: including benson14_varea is required to
                        % be able to obtain benson14_eccen, angle and sigma
specs.plotmesh      = 'left'; % left, right, both, or none
specs.plotelecs     = 'no'; % yes or no

out = electrode_to_nearest_node(specs);

%% to change one of the output figs to an 'empty' brain:
cmap = [1 1 1] * 0.7;colormap(cmap);colorbar off

%% to change view to lateral:
view(-90,0);

%% NEED TO FIX:UMCU

% UMCU patient electrode locations seem incorrect??

% specs = [];
% specs.pID      = 'beilen';
% specs.elecFile = '/Volumes/server/Projects/BAIR/Data/BIDS/visual/sub-umcubeilen/ses-umcuecogday03/ieeg/sub-beilen_ses-umcuecogday03_acq-clinical_electrodes.tsv';
% specs.fsDir    = '/Volumes/server/Freesurfer_subjects/umcu_beilen';
% specs.thresh   = []; 
% specs.patientPool = 'BAIR';
