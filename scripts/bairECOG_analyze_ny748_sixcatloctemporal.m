tbUse temporalECoG;

%% get data
projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual'; 
subject           = 'som748';
inputFolder       = 'ECoGBroadband';
description       = 'broadband';
session           = 'nyuecog01';
task              = [];
runnums           = [];

% set parameters for epoching
epoch_t           = [-0.2 1]; % epoch time window
chan_names        = {'GB'};   % which channels to include
stim_names        = {'ONEPULSE-1', 'ONEPULSE-2', 'ONEPULSE-3', 'ONEPULSE-4', 'ONEPULSE-5', 'ONEPULSE-6', ...
                     'TWOPULSE-1', 'TWOPULSE-2', 'TWOPULSE-3', 'TWOPULSE-4', 'TWOPULSE-5', 'TWOPULSE-6'};
                    % how to group the stimuli

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

% to add: 
% ecog_selectEpochs (needs to be improved)
opts = [];
opts.stimnames = stim_names;
opts.stimnames = stim_names;
opts.elec_selection_method = 'splithalf';
opts.elec_splithalf_thresh = 0.2;
[epochs, channels, chan_idx, R2, epochs_split] = ecog_selectElectrodes(epochs, channels, events, t, opts);

[epochs]    = ecog_averageEpochs(epochs, events, stim_names); 


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

stimcond = [1 1 1 1 1 1 2 2 2 2 2 2];
[results] = tde_evaluateModelFit(epochs, modelfun, params, pred, stimcond);

saveDir = '/Volumes/server/Projects/BAIR/Analyses/visual/sub-som748/ses-nyuecog01/figures';


tde_plotDataAndFits(results, epochs, channels, stim_ts, stim_info, t, {'ONEPULSE', 'TWOPULSE'}, saveDir);


tde_plotParams(results, channels, saveDir);

tde_plotDerivedPredictions(results,channels,2,1, saveDir);
