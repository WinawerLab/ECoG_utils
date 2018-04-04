
%tbUse('ECoG_utils');

%%%% DEFINE PATHS TO DATA, ANALYSIS DIRECTORIES %%%%

dataDir = '/Volumes/server/Projects/BAIR/Data/';
anaDir = '/Volumes/server/Projects/BAIR/Analyses/';
projectName = 'visual';

sub_label = 'beilen'; 
ses_label = 'day03';

allTasks = {'bairhrfpattern',...
            'bairprf',...
            'bairspatialobject',...
            'bairspatialpattern', ...
            'bairtemporalpattern'
            };
allRuns  = {'00', '01'};
numberofRunsPerTask = [1 2 2 2 2];

% Got this information from Gio (should perhaps be added to JSON files?)
n_chan = 133;
s_freq = 2048;

% NOTES:

%% Import ECoG data and Events

% Put all data in one struct
data = struct();
data.sub_label  = sub_label;
data.ses_label  = ses_label;
data.task_label = allTasks;
data.run_label  = allRuns;
data.nruns      = numberofRunsPerTask;

for cTask = 1:length(allTasks)
    
    task_label = allTasks{cTask};
    numberofRuns = numberofRunsPerTask(cTask);
    
    for cRun = 1:numberofRuns
        run_label = allRuns{cRun};
        
        dataToLoad = fullfile(dataDir, projectName, ['sub-' sub_label],['ses-' ses_label], 'ieeg', ...
            ['sub-' sub_label '_' 'ses-' ses_label '_' 'task-' task_label '_' 'run-' run_label]);

        dataName     =  [dataToLoad '_acq-clinical_ieeg.bin'];
        jsonName     =  [dataToLoad '_acq-clinical_ieeg.json'];
        eventsName   =  [dataToLoad '_acq-clinical_events.tsv'];       
        channelsName =  [dataToLoad '_acq-clinical_channels.tsv'];
        
        % Read in JSON file to get  duration for this run
        json = jsonread(jsonName);
        
        % Read in ECoG data
        duration = json.RecordingDuration; disp(duration);
        mem = memmapfile(dataName,  'Format', {'double', [n_chan duration * s_freq], 'data'});

        % Read in event info
        events = readtable(eventsName, 'FileType', 'text');
        
        % Read in channel info
        chans = readtable(channelsName, 'FileType', 'text');
        
        % Put all data in one struct
        data.dataset{cTask,cRun} = dataName;
        data.trial{cTask,cRun}   = mem.Data.data;
        data.label{cTask,cRun}   = chans;
        data.events{cTask,cRun}  = events;
        data.time{cTask,cRun}    = 0:1/s_freq:duration-(1/s_freq);
    end
end

data.fsample = s_freq;

%% PLOT: Powerspectra

% Pick a task and run (PROBABLY BETTER TO DO THIS ON CONCATENATED DATA)

task_to_plot = 'bairspatialpattern';
runidx = 1;

taskidx = find(strcmp(task_to_plot,data.task_label));

% Determine which channels are EEG channels
ecog_chans = find(contains(data.label{taskidx, runidx}.type, 'ECOG'));

% Two types of power spectra plots:
figure; spectopo(data.trial{taskidx,runidx}(ecog_chans,:),size(data.trial{taskidx,runidx},2),data.fsample, 'limits', [0 1024]);
[pxx,f] = pwelch(data.trial{taskidx,runidx}(ecog_chans,:)',data.fsample,0,data.fsample,data.fsample);
figure,plot(f,log10(pxx));

% Add onsets to indicate relevant data segments
onsets = data.events{taskidx,runidx}.onset;
    
% Plot variance per channel
figure; plot(1:length(ecog_chans), var(data.trial{taskidx,runidx}(ecog_chans,:),0,2), '-o');
set(gca, 'XTick', 1:length(ecog_chans), 'XTickLabel', data.label{taskidx,runidx}.name);

% Identify bad channels MANUALLY
bad_channels = [4 9 11 12 13 20]; % Manually selected based on inspection of spectra/time courses/variance plot
good_channels = setdiff(1:length(ecog_chans),bad_channels);

% % To inspect data on a channel by channel basis:
% for ii = 1:length(ecog_chans)
%     figure;hold on;
%     plot(data.time{taskidx,runidx}, data.trial{taskidx,runidx}(ii,:)); 
%     title([num2str(ii) ' ' data.label{taskidx,runidx}.name{ii}]);   
%     plot(onsets, zeros(length(onsets),1), 'ok')
%     waitforbuttonpress; close
% end

%% Visualize the timeseries of some channels
figure;
subplot(2,1,1); hold on
plot(data.time{taskidx,runidx},data.trial{taskidx,runidx}(bad_channels(1),:),'r')
title('first bad channel')
subplot(2,1,2); hold on
plot(data.time{taskidx,runidx},data.trial{taskidx,runidx}(good_channels(1),:),'k')
title('first good channel')

% Check the powerplot of all the good channels, no leftover outliers?
figure,plot(f,log10(pxx(:,good_channels)))

% Check the timeseries of all the good channels, no noisy moments?
figure,hold on;
plot(data.time{taskidx,runidx},data.trial{taskidx,runidx}(good_channels,:))
plot(onsets, zeros(length(onsets),1), 'ok')

%%  Do a common average reference using the good channels

plotCar = 'yes';

for cTask = 1:length(allTasks) 
    for cRun = 1:data.nruns(cTask)
        
        data.signal{cTask,cRun} = ecog_CarRegress(data.trial{cTask,cRun}, good_channels);
        
        switch plotCar
            case 'yes'
                figure;
                subplot(1,2,1); hold on
                channel_plot = good_channels(1);
                plot(data.time{cTask,cRun},data.trial{cTask,cRun}(channel_plot,:),'k')
                plot(data.time{cTask,cRun},data.signal{cTask,cRun}(channel_plot,:),'g')
                legend({'before CAR','after CAR'})
                title([data.task_label{cTask} ' ' data.run_label{cRun}]);
                
                subplot(1,2,2),hold on
                [pxx2,f] = pwelch(data.signal{cTask,cRun}(ecog_chans,:)',data.fsample,0,data.fsample,data.fsample);
                plot(f,pxx(:,channel_plot),'k')
                plot(f,pxx2(:,channel_plot),'g'); 
                set(gca, 'YScale', 'log')
        end
    end
end

%% Look at the effect of CAR



%% Save preprocessed data for further analysis

dataToSave = fullfile(anaDir, projectName, ['sub-' sub_label], ...
    ['sub-' sub_label '_' 'ses-' ses_label '_' 'task-' task_label '_preproc']);

save(dataToSave, 'data', 'signal', 'events');

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
plot(events.onset, zeros(length(events.onset),1),'k.','MarkerSize', 25, 'LineStyle','none');
legend(data.label(el));
xlabel('time (s)');
ylabel('broadband power (60-120 Hz)');
%ylabel('evoked');
l = line([data.time{1}(1) max(data.time{1})], [0 0],'LineStyle', '-', 'Color', 'k');
l.Annotation.LegendInformation.IconDisplayStyle = 'off';

%% Epoch data (all channels, all tasks)

% set epoch length (in seconds)
prestim = 0.5; 
poststim = 1.5;

% determine how many samples to go back and forth to extract epoch
onset_pre = round(prestim/(1/data.fsample));
onset_post = round(poststim/(1/data.fsample));

trials = struct();
    
[~,onsetsInx] = intersect(data.time{1},events.onset);

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
trials.events = events;
trials.fsample = data.fsample;

%% Save filtered and epoched data for further analysis

dataToSave = fullfile(anaDir, projectName, ['sub-' sub_label], ...
    ['sub-' sub_label '_' 'ses-' ses_label '_' 'task-' task_label '_epoched']);

save(dataToSave, 'trials'); % add index of interesting electrodes?




