% close all;

tbUse('ECoG_utils');

% check paths
if isempty(which('mrVista')), tbUse('vistasoft'); end


%% run script for a specific patient (with plotmesh on) to look at electrodes on brain

% Inputting entire paths for now until folder structure is settled

% specs.pID      = 'ny648';
% specs.elecFile = '/Volumes/server/Projects/BAIR/Data/visual/sub-ny648/ses-NYUECOG01/ieeg/sub-ny648_ses-NYUECOG01_electrodes.tsv';
% specs.fsDir    = '/Volumes/server/Freesurfer_subjects/som648';
% specs.thresh   = []; % default is 20 mm

specs.pID      = 'beilen';
specs.elecFile = '/Volumes/server/Projects/BAIR/Data/visual/sub-beilen/ses-day03/sub-beilen_acq-clinicalprojectedregions_electrodes.tsv';
specs.fsDir    = '/Volumes/server/Freesurfer_subjects/umcu_beilen';
specs.thresh   = []; % default is 20 mm

plotmesh  = 'yes'; % plot meshes with atlases for each subject: yes or no
plotlabel = 'yes'; % plot electrode labels on mesh: yes or no

retinotopicElecs = electrode_to_nearest_node(specs, plotmesh, plotlabel);




