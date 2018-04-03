
tbUse('ECoG_utils');

%%%% DEFINE PATHS TO DATA, ANALYSIS DIRECTORIES %%%%

dataDir = '/Volumes/server/Projects/BAIR/Data/';
anaDir = '/Volumes/server/Projects/BAIR/Analysis/';
projectName = 'visual';

sub_label = 'ny648'; 
ses_label = 'NYUECOG01';
task_label = 'soc';
run_label = '01';

% NOTES:
% - channel with triggers = 138 (DC10_REF)
% - note we also have HRF data for this subject; BIDS events file for these
% trials needs to be created

%% Import ECoG data and Events

dataToLoad = fullfile(dataDir, projectName, ['sub-' sub_label],['ses-' ses_label], 'ieeg', ...
    ['sub-' sub_label '_' 'ses-' ses_label '_' 'task-' task_label '_' 'run-' run_label]);

dataName =      [dataToLoad '_ieeg.edf'];
eventsName =    [dataToLoad '_events.tsv'];

% Read in ECoG data
cfg            = [];
cfg.dataset    = dataName;
cfg.continuous = 'yes';
cfg.channel    = 'all';
data           = ft_preprocessing(cfg);

% Read in events data
events         = readtable(eventsName, 'FileType', 'text');

%% Look at powerpectrum of channels to see which ones look bad

% Look at data.label to determine which channels are EEG channels: 1:126
chans = contains(data.label, 'EEG');

% This is one way to calculate the power spectrum and make a plot, the nice
% thing is that when you click on the lines it returns the channel numbers
figure;
spectopo(data.trial{1}(chans,:),size(data.trial{1},2),data.fsample);

% This is another way to look at the power spectum and make a plot. This
% returns the frequency (f) and power (pxx) and you can check for outliers.
[pxx,f] = pwelch(data.trial{1}(chans,:)',data.fsample,0,data.fsample,data.fsample);
figure,plot(f,log10(pxx))

% Identify bad channels MANUALLY
bad_channels = [81 82 100 102 126]; % Manually selected based on inspection of spectra and time courses
good_channels = setdiff(find(chans),bad_channels);

%% Visualize the timeseries of some channels

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

% TO ADD: plot events on top; check of onsets??
% TO ADD: extract electrode positions using runE2N, make list of visual electrodes??

%%  Do a common average reference using the good channels

signal = data.trial{1};
[signal] = ecog_CarRegress(signal, good_channels);

%% Look at the effect of CAR

figure;
subplot(1,2,1); hold on
channel_plot = good_channels(1);
plot(data.time{1},data.trial{1}(channel_plot,:),'k')
plot(data.time{1},signal(channel_plot,:),'g')
legend({'before CAR','after CAR'})

subplot(1,2,2),hold on
[pxx2,f] = pwelch(signal',data.fsample,0,data.fsample,data.fsample);
plot(f,pxx(:,channel_plot),'k')
plot(f,pxx2(:,channel_plot),'g'); 
set(gca, 'YScale', 'log')

%% Save preprocessed data for further analysis

dataToSave = fullfile(anaDir, projectName, ['sub-' sub_label], ...
    ['sub-' sub_label '_' 'ses-' ses_label '_' 'task-' task_label '_preproc']);

save(dataToSave, 'data', 'signal', 'events');

%% NEXT: nyuECOG_filter




