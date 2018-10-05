tbUse('ECoG_utils');

% SCRIPT DESCRIPTION %
% This script preprocesses the BIDS-formatted BAIR ECoG data into an evoked
% and a broadband response for each individual stimulus. It is meant be run
% cell-by-cell because some manual inputs are required (for channel
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
% Remarks: 
%
% Broadband computation: now done across entire concatenated time course,
% but could also be done for segments (but may introduce edge artifacts)
%
% Include baseline correction here?
% Add electrode plotting here?
% What can we take as 'no stimulation' for spectra plots? -> PRF = blanks,
% how about HRF/SOC?

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

%% [2] INSPECTION: Look at powerpectrum of channels to see which ones look bad

% Look at data.label to determine which channels are EEG channels: 1:126
chan_select = contains(channels.type, 'seeg');

% This is one way to calculate the power spectrum and make a plot, the nice
% thing is that when you click on the lines it returns the channel numbers
figure;
spectopo(data.trial{1}(chan_select,:),size(data.trial{1},2),data.fsample);

% This is another way to look at the power spectum and make a plot. This
% returns the frequency (f) and power (pxx) and you can check for outliers.
[pxx,f] = pwelch(data.trial{1}(chan_select,:)',data.fsample,0,data.fsample,data.fsample);
figure,plot(f,log10(pxx));
xlabel('Frequency (Hz)'); ylabel('Log power');

% Identify channels whose average logpower across frequencies is 2 sd
% above/below the mean across channels:
mn = mean(log10(pxx),1);
sd = std(mn,0,2);
disp('Two SD above the mean:'); 
disp(find(mn>(mean(mn,2)+2*sd)));
disp('Two SD below the mean:');
disp(find(mn<(mean(mn,2)-2*sd)));
% !SD should not be used to automatically identify bad channels, because
% channels with strong activation can have higher power on average!
% Instead, pay extra attention to those channels when going through them:

% Look at the entire time course for each channel:
for cChan = 1:size(data.trial{1},1); figure;plot(data.trial{1}(cChan,:)); title([num2str(cChan) ': ' data.label{cChan}]);waitforbuttonpress;close; end

%% [2] INSPECTION: Identify bad channels here:

bad_channels = [100 102 126]; % Manually selected based on inspection of spectra and time courses
good_channels = setdiff(find(chan_select),bad_channels);

%% [2] INSPECTION: Visualize the timeseries of some channels

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

%% [3] PREPROCESS: Do a common average reference using the good channels

signal = data.trial{1};
[signal] = ecog_CarRegress(signal, good_channels);

%% [3] INSPECTION: Look at the effect of CAR

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

%% [4] PREPROCESS: Compute time-varying broadband 

% Define frequency bands to filter
bands = [[70 90]; [90 110]; [130 150]; [150 170]];

% Define broadband extraction method (see help extractBroadband_nyu)
bbmethod = 7;

% Extract broadband
broadband = extractBroadband_nyu(signal', data.fsample, bbmethod, bands);
hband_sig = broadband';

%% [4] INSPECTION: Plot filtered time course and broadband for channel(s) of interest

% electrodes with coverage in visual areas: 
eltomatch = {'MO_01'}; %, 'MO_02', 'MO_03', 'MO_04'};

% NY648 electrodes with coverage in IPS:
%eltomatch = {'DLPA_07', 'DLPA_08'};
%eltomatch = {'G_03', 'G_04', 'G_09', 'G_10', 'G_11', 'G_12', 'G_17', 'G_18'};
%eltomatch = {'G_08', 'G_07'};

% find matching electrode numbers
el = ecog_matchchannels(eltomatch, data);

% smooth & plot
figure; hold on
for ii = 1:length(eltomatch)
    plot(data.time{1},smooth(hband_sig(el(ii),:),128), 'LineWidth', 2);
    %plot(data.time{1},smooth(signal(el(ii),:),128), 'LineWidth', 2);
end

% and plot event onsets on top
plot(eventsTable.onset, zeros(length(eventsTable.onset),1),'k.','MarkerSize', 25, 'LineStyle','none');
legend(data.label(el));
xlabel('time (s)');
ylabel('broadband power (60-120 Hz)');
%ylabel('evoked');
l = line([data.time{1}(1) max(data.time{1})], [0 0],'LineStyle', '-', 'Color', 'k');
l.Annotation.LegendInformation.IconDisplayStyle = 'off';

%% [5] Segment the data (all channels, all tasks)

% TO DO include NO STIM segments --> use HRF runs e.g. 1.5 to 0.5 seconds prestim

% Segmentation
%epochLength = 1.5; % in seconds
%baselineLength = -0.5; % in seconds

%eventsInSamples = events.onset*hdr.Fs;
    %cfg.trl{iRun} = [eventsInSamples eventsInSamples+(epochLength*hdr.Fs) repmat(baselineLength*hdr.Fs, size(eventsInSamples))];

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

