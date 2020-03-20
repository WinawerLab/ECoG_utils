
clear specs;

projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual'; 
subject           = 'som748';
session           = 'nyuecog04';
task              = 'sixcatloctemporal';
%session           = 'nyuecog01';
%task              = 'temporalpattern';
runnums           = '01';
%specs.stim_names  = {'BODIES', 'BUILDINGS', 'FACES', 'OBJECTS', 'SCENES', 'SCRAMBLED'};

specs.stim_names  = {'TWOPULSE-1', 'TWOPULSE-2', 'TWOPULSE-3', 'TWOPULSE-4', 'TWOPULSE-5', 'TWOPULSE-6'};
%specs.stim_names  = {'TWOPULSE-1-0', 'TWOPULSE-2-0', 'TWOPULSE-3-0', 'TWOPULSE-4-0', 'TWOPULSE-5-0', 'TWOPULSE-6-0'};
%specs.stim_names  = {'TWOPULSE-1-1', 'TWOPULSE-2-1', 'TWOPULSE-3-1', 'TWOPULSE-4-1', 'TWOPULSE-5-1', 'TWOPULSE-6-1'};
%specs.stim_names  = {'ONEPULSE-1', 'ONEPULSE-2', 'ONEPULSE-3', 'ONEPULSE-4', 'ONEPULSE-5', 'ONEPULSE-6'};

%specs.stim_names  = {'TWOPULSE-1-FACES', 'TWOPULSE-2-FACES', 'TWOPULSE-3-FACES', 'TWOPULSE-4-FACES', 'TWOPULSE-5-FACES', 'TWOPULSE-6-FACES'};
%specs.stim_names  = {'ONEPULSE-1-FACES', 'ONEPULSE-2-FACES', 'ONEPULSE-3-FACES', 'ONEPULSE-4-FACES', 'ONEPULSE-5-FACES', 'ONEPULSE-6-FACES'};
specs.chan_names  = {'GB018', 'GB034'};
specs.plot_ylim   = [-1 6];
specs.plot_cmap   = 'copper';
specs.plot_includelegend = 0;

bidsEcogPlotTrials(projectDir, subject, session, task, runnums, [], [], specs, 1);
%bidsEcogPlotTrials(projectDir, subject, session, task);

%%

tbUse temporalECoG;


%% get data
projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual'; 
subject           = 'som748';
inputFolder       = 'ECoGBroadband';
description       = 'broadband';
session           = 'nyuecog04';
task              = [];
runnums           = '01';

% set parameters for epoching
epoch_t           = [-0.2 1]; % epoch time window
chan_names        = {'GB'};   % which channels to include
stim_names        = {'ONEPULSE-1', 'ONEPULSE-2', 'ONEPULSE-3', 'ONEPULSE-4', 'ONEPULSE-5', 'ONEPULSE-6', ...
                     'TWOPULSE-1', 'TWOPULSE-2', 'TWOPULSE-3', 'TWOPULSE-4', 'TWOPULSE-5', 'TWOPULSE-6'};
                    % how to group the stimuli
stim_names        = {'BODIES', 'BUILDINGS', 'FACES', 'OBJECTS', 'SCENES', 'SCRAMBLED'};

% load data
dataPath = fullfile(projectDir, 'derivatives', inputFolder);
[data, channels, events] = bidsEcogGetPreprocData(dataPath, subject, session, task, runnums, description);

% select channels 
chan_idx    = contains(channels.name, chan_names); 
data        = data(chan_idx,:);
channels    = channels(chan_idx,:);

% select stimuli
stim_idx    = contains(events.trial_name, stim_names);
events      = events(stim_idx,:);

% make epochs
[epochs, t] = ecog_makeEpochs(data, events.onset, epoch_t, channels.sampling_frequency(1));  
[epochs]    = ecog_normalizeEpochs(epochs, t, [min(epoch_t) 0], 'percentsignalchange');
[epochs]    = ecog_averageEpochs(epochs, events, stim_names); 

% to add: 
% ecog_selectEpochs (needs to be improved)
% ecog_selectElectrodes (needs to be written based on code in tde_selectData)

% get visual area match per electrode
opt = [];
opt.pID           = subject; 
opt.plotmesh      = 'none';
visualelectrodes  = electrode_to_nearest_node(opt);

% add visual areas to channel table
[channels] = bair_addVisualAtlasNamesToChannelTable(channels,visualelectrodes);

%% fit model

[stim_ts, stim_info] = tde_generateStimulusTimecourses(stim_names,t);

modelfuns = tde_modelTypes();
modelfun  = modelfuns(1); 

options = struct();
options.display  = 'off';

[params, pred] = tde_fitModel(modelfun, stim_ts, epochs, 512, options);

%% visualize model fits

[results] = tde_evaluateModelFit(epochs, modelfun, params, pred);

saveDir = '/Volumes/server/Projects/BAIR/Analyses/visual/sub-som748/ses-nyuecog04/figures';

tde_plotDataAndFits(results, epochs, channels, stim_ts, stim_info, t, {'ONEPULSE', 'TWOPULSE'}, saveDir);

tde_plotFittedAndDerivedParams(results, channels, saveDir);
