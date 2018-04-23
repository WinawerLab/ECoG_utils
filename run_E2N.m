% close all;
tbUse('ECoG_utils');

%% run script for a specific patient (with plotmesh on) to look at electrodes on brain

% Inputting entire paths for now until folder structure is settled

specs.pID      = 'ny648';
specs.elecFile = '/Volumes/server/Projects/BAIR/Data/BIDS/visual/sub-ny648/ses-NYUECOG01/ieeg/sub-ny648_ses-NYUECOG01_electrodes.tsv';
specs.fsDir    = '/Volumes/server/Freesurfer_subjects/som648';
specs.thresh   = []; 
specs.patientPool = 'BAIR';

% specs.pID      = 'beilen';
% specs.elecFile = '/Volumes/server/Projects/BAIR/Data/BIDS/visual/sub-beilen/ses-day03/sub-beilen_acq-clinicalprojectedregions_electrodes.tsv';
% specs.fsDir    = '/Volumes/server/Freesurfer_subjects/umcu_beilen';
% specs.thresh   = []; 

% specs.pID      = '668';
% specs.elecFile = [];
% specs.fsDir    = '/Volumes/server/Freesurfer_subjects/som668';
% specs.thresh   = []; 
% specs.patientPool = 'SOM';

specs.plotmesh  = 'yes'; % plot meshes with atlases for each subject: yes or no
specs.plotlabel = 'yes'; % plot electrode labels on mesh: yes or no

out = electrode_to_nearest_node(specs);

%% PATH TO NOAH COLORMAPS:
% '/Volumes/server/Projects/HCP/analysis/images'


% alignment to fs average
% find nearest node ID in pial surface of indiv.subj
% find coordinate of this node ID in indiv.subj lh.sphere.reg
% find nearest node ID with closest coordinate in fsaverage lh.sphere
% that node ID gives you coordinates on lh.pial of fsaverage

