projectDir = '/Volumes/server/Projects/BAIR/Data/BIDS/tactile'; 
subject = 'umcudrouwen';

%% common average reference
bidsEcogRereference(projectDir, subject);

%% broadband 
%bands             = [[60 70]; [70 80]; [80 90]; [110 120]; [120 130]; [130 140]; [160 170]; [170 180]; [180 190]];
bands             = [[60 70]; [70 80]; [80 90]; [120 130]; [130 140]; [160 170]; [170 180]; [180 190]];
bidsEcogBroadband(projectDir, subject, [], [], [], bands);

%% plot
clear specs;

session           = 'umcuiemu01';
task              = 'vtstemporalpattern';
%inputFolder       = 'ECoGBroadband';
%description       = 'broadband';

inputFolder       = 'ECoGCAR';
description       = 'reref';

%specs.stim_names  = {'thumb', 'index', 'middle', 'ring', 'little'};
specs.epoch_t     = [-0.5 0.5];
specs.base_t      = [-0.5 -0.1];
specs.plot_type   = 'singletrial';
specs.stim_names  = {'ONEPULSE-1', 'ONEPULSE-2', 'ONEPULSE-3', 'ONEPULSE-4', 'ONEPULSE-5', 'ONEPULSE-6'};
%specs.stim_names  = {'TWOPULSE-1', 'TWOPULSE-2', 'TWOPULSE-3', 'TWOPULSE-4', 'TWOPULSE-5', 'TWOPULSE-6'};

%specs.chan_names  = {'C03', 'C04', 'C11', 'C12'};
specs.chan_names  = {'C14', 'C15', 'C17', 'C18'};

% specs.chan_names  = {'C'};
% specs.subplotdims = [4 8];
% specs.subplotidx  = [32:-1:1];

%specs.plot_ylim   = [-1 6];

savePlot = 1;

bidsEcogPlotTrials(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot);


%%

clear specs;

session           = 'umcuiemu01';
task              = 'vtstemporalpattern';

%specs.epoch_t     = [-0.5 0.5];
%specs.plot_type   = 'singletrial';
specs.stim_names  = {'ONEPULSE-1', 'ONEPULSE-2', 'ONEPULSE-3', 'ONEPULSE-4', 'ONEPULSE-5', 'ONEPULSE-6'};
%specs.stim_names  = {'TWOPULSE-1', 'TWOPULSE-2', 'TWOPULSE-3', 'TWOPULSE-4', 'TWOPULSE-5', 'TWOPULSE-6'};

%specs.chan_names  = {'C03', 'C04', 'C11', 'C12'};
%specs.chan_names  = {'C14', 'C15', 'C17', 'C18'};

specs.chan_names  = {'C'};
specs.subplotdims = [4 8];
specs.subplotidx  = [32:-1:1];

%specs.plot_ylim   = [-1 6];

savePlot = 1;

bidsEcogPlotSpectra(projectDir, subject, session, task, [], [], [], specs, savePlot);


