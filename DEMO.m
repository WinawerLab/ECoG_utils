% START: Convert Raw data received from SoM with bids-formatted data using
edit bairECOG_convertToBIDS_TEMPLATE.m
% Save script with patient- and session-specific extension patient (e.g.,
% bairECOG_convertToBIDS_ny748_ses1.m) and update specs and manual sections
% with patient relevant info, run it to convert data to BIDS.

% Once data is in BIDS, the functions below can be run to process the data.

%% derivatives 1: rereference 
projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual'; 
subject           = 'som748';
session           = 'nyuecog01';
task              = 'temporalpattern';

bidsEcogRereference(projectDir, subject, session, task)

%% derivatives 2: broadband computation
bidsEcogBroadband(projectDir, subject, session, task)

%% derivatives 3: match electrodes to atlases

projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual'; 
subject           = 'chaam';
atlasName         = {'benson14_varea', 'wang15_fplbl'};

electrode_table   = bidsEcogMatchElectrodesToAtlas(projectDir, subject, [], atlasName);
% TO DO make it write out electrode_table to derivatives folder

%% plots 1: visualizing electrodes on mesh (OLD)
opt = [];

opt.pID           = 'som748'; % patient ID 
opt.atlasNames    =  {'benson14_varea', 'wang15_mplbl'};
opt.plotmesh      = 'left';
opt.plotelecs     = 'yes';
opt.plotlabel     = 'yes';
out = electrode_to_nearest_node(opt,1);

% TO DO replace with bidsEcogPlotElectrodesOnMesh.m

%% plots 2: time-varying responses

inputFolder       = 'ECoGBroadband';
description       = 'broadband';
runnums           = '01';

specs = [];
specs.epoch_t     = [-0.2 1.5];
specs.chan_names  = {'GB'};
specs.plot_type   = 'averageSE';
specs.stim_names  = {'TWOPULSE-6'};
specs.plot_ylim   = [-1 5];
    
bidsEcogPlotTrials(projectDir, subject, sessions, tasks, runnums, ...
        inputFolder, description, specs, savePlot) 
    
%% plots 3: response spectra

bidsEcogPlotSpectra(projectDir, subject, sessions, tasks, runnums, ...
    inputFolder, description, specs, savePlot)

