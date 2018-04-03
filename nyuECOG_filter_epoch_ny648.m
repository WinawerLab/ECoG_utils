tbUse('ECoG_utils');

%%%% DEFINE PATHS TO DATA, ANALYSIS DIRECTORIES %%%%

anaDir = '/Volumes/server/Projects/BAIR/Analysis/';
projectName = 'visual';

sub_label = 'ny648'; 
ses_label = 'NYUECOG01';
task_label = 'soc';

%% Load preprocessed data

dataToLoad = fullfile(anaDir, projectName, ['sub-' sub_label], ...
    ['sub-' sub_label '_' 'ses-' ses_label '_' 'task-' task_label '_preproc']);

load(dataToLoad);

%% Filter to extract broadband responses

% Define frequency bands to filter
bands = [[70 90]; [90 110]; [130 150]; [150 170]];
% bands = [[70 80]; [80 90]; [90 100]; [100 110]; [130 140]; [140 150]; [150 160]; [160 170]];

% Define broadband extraction method (see help extractBroadband_nyu)
bbmethod = 5;

% Extract broadband
broadband = extractBroadband_nyu(signal', data.fsample, bbmethod, bands);
hband_sig = broadband';

%% Plot filtered time course for channel(s) of interest
% load(stimFile); replace with events

% electrodes with coverage in visual areas: 
eltomatch = {'MO_01'}; %, 'MO_02', 'MO_03', 'MO_04'};

% NY648 electrodes with coverage in IPS:
%eltomatch = {'DLPA_07', 'DLPA_08'};
%eltomatch = {'G_03', 'G_04', 'G_09', 'G_10', 'G_11', 'G_12', 'G_17', 'G_18'};
%eltomatch = {'G_08', 'G_07'};

% find matching electrode numbers
el = ecog_matchchannels(eltomatch, data);

load(stimFile);
tasklabels = fields(onsets);

for tl = 1:length(tasklabels)
    tasklabel = tasklabels{tl};
    
    % smooth & plot
    figure('Name', tasklabel), hold on,
    for ii = 1:length(eltomatch)
        plot(data.time{1},smooth(hband_sig(el(ii),:),128), 'LineWidth', 2);
        %plot(data.time{1},smooth(signal(el(ii),:),128), 'LineWidth', 2);
    end

    % and plot event onsets on top
    plot(onsets.(tasklabel), zeros(length(onsets.(tasklabel)),1),'k.','MarkerSize', 25, 'LineStyle','none');
    legend(data.label(el));
    xlabel('time (s)');
    ylabel('broadband power (60-120 Hz)');
    %ylabel('evoked');
    l = line([data.time{1}(1) max(data.time{1})], [0 0],'LineStyle', '-', 'Color', 'k');
    l.Annotation.LegendInformation.IconDisplayStyle = 'off';
end

%% Epoch data (all channels, all tasks)

% set epoch length (in seconds)
prestim = 0.5; 
poststim = 1.5;

% determine how many samples to go back and forth to extract epoch
onset_pre = round(prestim/(1/data.fsample));
onset_post = round(poststim/(1/data.fsample));

trials = struct();
for tl = 1:length(tasklabels)
    
    tasklabel = tasklabels{tl};
    [~,onsetsInx] = intersect(data.time{1},onsets.(tasklabel));

    % extract epochs
    for ii = 1:length(data.label)
        for jj = 1:length(onsetsInx)  
            % broadband
            trials.broadband.(tasklabel)(ii,:,jj) = smooth(hband_sig(ii,onsetsInx(jj)-onset_pre:onsetsInx(jj)+onset_post),8);
            % evoked
            trials.evoked.(tasklabel)(ii,:,jj) = smooth(signal(ii,onsetsInx(jj)-onset_pre:onsetsInx(jj)+onset_post),8);
        end
    end
end
disp('done');

trials.label = data.label;
trials.time = -prestim:(1/data.fsample):poststim;
trials.stimuli = stimuli;
trials.fsample = data.fsample;

%% Save filtered and epoched data for further analysis

dataToSave = fullfile(anaDir, projectName, ['sub-' sub_label], ...
    ['sub-' sub_label '_' 'ses-' ses_label '_' 'task-' task_label '_epoched']);

save(dataToSave, 'trials'); % add index of interesting electrodes?


