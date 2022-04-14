
% BIDS formatted iEEG data can be processed using the functions below:

%% select a BIDS project, subject, session and task
projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual_ecog_recoded'; %% NEEDS TO BE CHANGED TO POINT TO OPENNEURO 
subject           = 'p01';   
sessions          = 'umcuiemu01';
tasks             = 'temporalpattern';

%% derivatives 1: rereference 
bidsEcogRereference(projectDir, subject, sessions, tasks)

%% derivatives 2: broadband computation
bidsEcogBroadband(projectDir, subject, sessions, tasks)

%% derivatives 3: match electrodes to atlases
atlasName         = {'wang15_mplbl'};
electrode_table   = bidsEcogMatchElectrodesToAtlas(projectDir, subject, [], atlasName);

%% plots 1: visualizing electrodes on mesh 
specs = [];
specs.plotelecs     = 'yes';
specs.plotlabel     = 'yes';
bidsEcogPlotElectrodesOnMesh(projectDir, subject, [], atlasName, [], specs);

%% plots 2: time-varying responses
inputFolder       = 'ECoGBroadband';
description       = 'broadband';
runnums           = '01';

specs = [];
specs.epoch_t     = [-0.2 1.5];
specs.chan_names  = {'OT08'};
specs.plot_type   = 'averageSE';
specs.stim_names  = {'TWOPULSE-6'};
specs.plot_ylim   = [-1 5];
    
savePlot          = true;

bidsEcogPlotTrials(projectDir, subject, sessions, tasks, runnums, ...
        inputFolder, description, specs, savePlot) 
    
%% plots 3: response spectra

bidsEcogPlotSpectra(projectDir, subject, sessions, tasks, runnums, ...
    inputFolder, description, specs, savePlot)

