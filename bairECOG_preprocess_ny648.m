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
%
% TO DO: add prestimulus baseline correction (here?)
% TO DO: add electrode plotting (here?), add list of visual electrodes to
% data outputs
% TO DO: add 'no stimulation' segments -> PRF = blanks, HRF take -1.5:-0.5

%% [0] Define paths and dataset specs %%

% Dataset specs
projectName = 'visual';
sub_label   = 'som648'; 
ses_label   = 'nyuecog01';
task_labels = {'hrfpattern','soc'};

% Input paths specs
dataPth     = '/Volumes/server/Projects/BAIR/Data/BIDS/';

% Output paths specs
dataDir = fullfile(dataPth, projectName, sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label), 'ieeg');
stimDir = fullfile(dataPth, projectName, 'stimuli');
saveDir = fullfile(dataPth, projectName, 'derivatives', 'preprocessed', sprintf('sub-%s', sub_label));

%% [1] Read in ECoG data

[data, events, channels]  = ecog_readBIDSData(dataDir, sub_label, ses_label);

%% [1] Plot electrode locations for this patient to see if there are matches with visual regions

specs = [];
specs.pID           = sub_label; % patient ID number
specs.patientPool   = 'BAIR';

visualelectrodes = electrode_to_nearest_node(specs, dataDir);

%% [2] DATA INSPECTION: Look at powerpectrum of channels to see which ones look bad

% Look at data.label to determine which channels are EEG channels: 1:126
chan_select = contains(channels.type, 'seeg');

% Generate spectral plot; check command window output for outliers 
outliers = ecog_plotChannelSpectra(data, chan_select);

% NOTE: outliers (identified as channels with mean power that is more that
% two standard deviations above or below the average across chanensl should
% not be used to automatically identify bad channels, because channels with
% strong activation can have higher power on average! Instead, pay extra
% attention to those channels when inspecting their time courses:

% Inspect the time course for suspicious / outlier channels
for channel_plot = 1:length(outliers); figure;plot(data.time{1}, data.trial{1}(outliers(channel_plot),:)); title([num2str(outliers(channel_plot)) ': ' data.label{outliers(channel_plot)}]); end

% Even better: look at the time course for each channel in succession:
for channel_plot = 1:size(data.trial{1},1); figure;plot(data.time{1}, data.trial{1}(channel_plot,:)); title([num2str(channel_plot) ': ' data.label{channel_plot}]);waitforbuttonpress;close; end

%% [2] Identify bad channels here:

bad_channels = [100 102 126]; % Manually selected based on inspection of spectra and time courses
good_channels = setdiff(find(chan_select),bad_channels);

%% [2] DATA INSPECTION: Visualize the timeseries of some channels

figure;
subplot(2,1,1); hold on
plot(data.time{1},data.trial{1}(bad_channels(1),:),'r')
title('first bad channel')
subplot(2,1,2); hold on
plot(data.time{1},data.trial{1}(good_channels(1),:),'k')
title('first good channel')

% Check the powerplot of all the good channels, no leftover outliers?
figure,plot(f,log10(pxx(:,good_channels)))

% Check the timeseries of all the good channels, no noisy moments?
figure,plot(data.time{1},data.trial{1}(good_channels,:))

%% [3] Do a common average reference using the good channels

signal = ecog_CarRegress(data.trial{1}, good_channels);

%% [3] DATA INSPECTION: Look at the effect of CAR

figure;
subplot(1,2,1); hold on
channel_plot = good_channels(1);
plot(data.time{1},data.trial{1}(channel_plot,:),'k')
plot(data.time{1},signal(channel_plot,:),'g')
legend({'before CAR','after CAR'})
xlabel('Time (s)'); ylabel('Voltage');
title(data.label(channel_plot));

subplot(1,2,2),hold on
[pxx2,f] = pwelch(signal',data.fsample,0,data.fsample,data.fsample);
plot(f,pxx(:,channel_plot),'k')
plot(f,pxx2(:,channel_plot),'g'); 
set(gca, 'YScale', 'log')
xlabel('Frequency (Hz)'); ylabel('Log power');
title(data.label(channel_plot));

%% [4] Compute time-varying broadband 

% Define frequency bands to filter
bands = [[70 90]; [90 110]; [130 150]; [150 170]];

% Define broadband extraction method (see help extractBroadband_nyu)
bbmethod = 7;

% Extract broadband
[broadband, methodstr] = extractBroadband_nyu(signal', data.fsample, bbmethod, bands);

%% [4] DATA INSPECTION: Plot filtered time course and broadband for channel(s) of interest

% Pick one or more of the electrodes with coverage in visual areas: 
disp(visualelectrodes.benson14_varea);

% e.g.
eltomatch = visualelectrodes.benson14_varea.elec_labels(9:12); %MO01-04

% Find matching channel number
el = ecog_matchChannels(eltomatch, data);

% Plot evoked and broadband
ecog_plotFullTimeCourse(timecourse,chanstoplot,events)

%% [5] Segment the data (all channels, all tasks)

% TO DO include NO STIM segments --> use HRF runs e.g. 1.5 to 0.5 seconds prestim
% TO DO no smoothing in preprocessing!

% Use fieldtrip trlfunction?
%epochLength = 1.5; % in seconds
%baselineLength = -0.5; % in seconds

eventsInSamples = events.onset*hdr.Fs;
cfg.trl = [eventsInSamples eventsInSamples+(epochLength*hdr.Fs) repmat(baselineLength*hdr.Fs, size(eventsInSamples))];

% set epoch length (in seconds)
prestim = 0.5; 
poststim = 1.5;

% determine how many samples to go back and forth to extract epoch
onset_pre = round(prestim/(1/data.fsample));
onset_post = round(poststim/(1/data.fsample));

trials = struct();
    
[~,onsetsInx] = intersect(data.time{1},eventsTable.onset);

% extract epochs
for ii = 1:length(data.label)
    for jj = 1:length(onsetsInx)  
        % broadband
        trials.broadband(ii,:,jj) = smooth(hband_sig(ii,onsetsInx(jj)-onset_pre:onsetsInx(jj)+onset_post),8);
        % evoked
        trials.evoked(ii,:,jj) = smooth(signal(ii,onsetsInx(jj)-onset_pre:onsetsInx(jj)+onset_post),8);
    end
end

disp('done');

trials.label = data.label;
trials.time = -prestim:(1/data.fsample):poststim;
trials.events = eventsTable;
trials.fsample = data.fsample;

%% Save filtered and epoched data for further analysis

dataToSave = fullfile(saveDir, ['sub-' sub_label '_preproc']);

save(dataToSave, 'data', 'trials'); % add index of interesting electrodes?

