tbUse('ECoG_utils');

% SCRIPT DESCRIPTION %
% This script preprocesses the BIDS-formatted BAIR ECoG data. It is meant
% be run cell-by-cell because some manual inputs are required (for channel
% selection) and inspection of data.
%
% This is the order of operations:
%
% [0] Path specifications: data should be in BIDS
% [1] Reading in data, events and channel info
% [2] Channel inspection, selection of bad channels (MANUAL)
% [3] Common average reference (regression based)
% [4] Broadband computation
% [5] Segmentation
% 
% OUTPUT:
% - data: struct with entire raw and broadband time courses. 
% - trials: struct with evoked and a broadband response for each stimulus. 
%
% Remarks: 
%
% Broadband computation: now done across entire concatenated time course,
% but could also be done for segments (but may introduce edge artifacts?)

%% [0] Define paths and dataset specs %%

% Dataset specs
projectName = 'visual';
sub_label   = 'som692'; 
ses_label   = 'nyuecog01';

% Input paths specs
dataPth     = '/Volumes/server/Projects/BAIR/Data/BIDS/';

% Output paths specs
dataDir = fullfile(dataPth, projectName, sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label), 'ieeg');
saveDir = fullfile(dataPth, projectName, 'derivatives', 'preprocessed', sprintf('sub-%s', sub_label));

if ~exist(saveDir, 'dir'); mkdir(saveDir);end

%% [1] Read in ECoG data

[ftdata, events, channels]  = ecog_readBIDSData(dataDir, sub_label, ses_label);

%% [1] Plot electrode locations for this patient to see if there are matches with visual regions

specs = [];
specs.pID           = sub_label; % patient ID number
specs.patientPool   = 'BAIR';

visualelectrodes = electrode_to_nearest_node(specs, dataDir);

%% [3] Do a common average reference using the good channels

good_channels = find(contains(channels.status, 'good'));
signal = ecog_CarRegress(ftdata.trial{1}, good_channels);

%% [3] DATA INSPECTION: Look at the effect of CAR

% figure;
% subplot(1,2,1); hold on
% channel_plot = good_channels(1);
% plot(ftdata.time{1},ftdata.trial{1}(channel_plot,:),'k')
% plot(ftdata.time{1},signal(channel_plot,:),'g')
% legend({'before CAR','after CAR'})
% xlabel('Time (s)'); ylabel('Voltage');
% title(ftdata.label(channel_plot));
% 
% subplot(1,2,2),hold on
% [pxx2,freqs] = pwelch(signal',ftdata.fsample,0,ftdata.fsample,ftdata.fsample);
% plot(freqs,pxx(:,channel_plot),'k')
% plot(freqs,pxx2(:,channel_plot),'g'); 
% set(gca, 'YScale', 'log')
% xlabel('Frequency (Hz)'); ylabel('Log power');
% title(ftdata.label(channel_plot));

%% [4] Compute time-varying broadband 

% Define frequency bands to filter
bands = [[70 90]; [90 110]; [130 150]; [150 170]];

% Define broadband extraction method (see help ecog_extractBroadband)
bbmethod = 7;

% Extract broadband
[broadband, methodstr] = ecog_extractBroadband(signal', ftdata.fsample, bbmethod, bands);

%% Create a data structure to save

data = struct();

data.hdr        = ftdata.hdr;
data.time       = ftdata.time{1};
data.raw        = ftdata.trial{1}; 
data.car_reref  = signal;
data.broadband  = broadband;
data.bb_bands   = bands;
data.bb_method  = methodstr;
data.events     = events;
data.channels   = channels;
data.viselec    = visualelectrodes;
data.cfg        = ftdata.cfg;

%% [4] DATA INSPECTION: 

%% Plot filtered time course and broadband for channel(s) of interest

% Pick one or more of the electrodes with coverage in visual areas: 
disp(data.viselec.benson14_varea);

% e.g.
eltomatch = visualelectrodes.benson14_varea.elec_labels(12:15); %MO01-04

% Find matching channel number
el = ecog_matchChannels(eltomatch, data);

% Plot voltage and broadband
%ecog_plotFullTimeCourse(data,'car_reref', el);
%ecog_plotFullTimeCourse(data,'broadband', el); % no smoothing
%ecog_plotFullTimeCourse(data,'broadband', el, data.hdr.Fs/2); % with smoothing

%% Check whether the onsets align with the triggers by plotting trigger channel itself

eltomatch = data.channels.name{find(contains(data.channels.type,'trig'))};
el = ecog_matchChannels(eltomatch, data);

ecog_plotFullTimeCourse(data,'raw', el); title('triggers');

%% [5] Segment the data (all channels, all tasks)

% Set epoch length (in seconds)
epoch = [-0.5 1.5];

% Epoch the stimulus trials
trials = ecog_epochData(data, epoch);

% Add 'no stimulation' baseline trials --> not necessary here because we have
% blanks in the PRF runs?

%mockdata = data;
%mockdata.events = mockdata.events(strmatch('HRFPATTERN',data.events.trial_name),:);

% Define an epoch in which there were no stimuli
%baseline_epoch = [-1.5 -0.5];

% Epoch the baseline 'trials'
%baseline_trials = ecog_epochData(mockdata, baseline_epoch);
%baseline_trials.events = [];

%% Save preprocessed and epoched data for further analysis

% disp('saving preprocessed data...');
% saveName = fullfile(saveDir, sprintf('sub-%s_ses-%s_preproc', sub_label, ses_label));
% save(saveName, 'data'); 

disp('saving epoched data...');
saveName = fullfile(saveDir, sprintf('sub-%s_ses-%s_epoched_onflips', sub_label, ses_label));
%save(saveName, 'trials', 'baseline_trials'); 
save(saveName, 'trials'); 
disp('done');