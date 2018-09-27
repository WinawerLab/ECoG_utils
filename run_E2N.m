% close all;
tbUse('ECoG_utils');

%% run script for a specific patient (with plotmesh on) to look at electrodes on brain

specs = [];

specs.pID           = '668'; % patient ID number
specs.atlasNames    = {'wang2015_atlas', 'benson14_varea', 'benson14_eccen'}; % default is all of these maps: 
                                                            % {'wang2015_atlas',
                                                            % 'benson14_varea',
                                                            % 'benson14_eccen',
                                                            % 'benson14_angle',
                                                            % 'benson14_sigma',
                                                            % 'template_areas'};
specs.plotmesh      = 'yes'; % plot meshes with atlases for each subject: yes or no
specs.plotlabel     = 'yes'; % plot electrode labels on mesh: yes or no
specs.plotcolorbar  = 'yes'; % plot colorbar

out = electrode_to_nearest_node(specs);

%% NEED TO FIX:UMCU

% UMCU patient electrode locations seem incorrect??

% specs = [];
% specs.pID      = 'beilen';
% specs.elecFile = '/Volumes/server/Projects/BAIR/Data/BIDS/visual/sub-umcubeilen/ses-umcuecogday03/ieeg/sub-beilen_ses-umcuecogday03_acq-clinical_electrodes.tsv';
% specs.fsDir    = '/Volumes/server/Freesurfer_subjects/umcu_beilen';
% specs.thresh   = []; 
% specs.patientPool = 'BAIR';
