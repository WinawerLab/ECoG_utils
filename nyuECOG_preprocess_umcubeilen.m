
tbUse('ECoG_utils');

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

%% PLOT: Visualize the timeseries of some channels

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

%%  Do a common average reference using the good channels (but leave other channels in the data)

plotCar = 'yes';

for cTask = 1:length(allTasks) 
    for cRun = 1:data.nruns(cTask)
        
        data.signal{cTask,cRun} = ecog_CarRegress(data.trial{cTask,cRun}, good_channels);
        
        % Look at the effect of CAR
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

%% Filter to extract broadband responses

% Define frequency bands to filter: remember LINE NOISE at 50Hz
bands = [[60 70]; [70 80]; [80 90]; [110 120]];
% bands = [[60 70]; [80 90]; [110 120]; [120 130]; [130 140]; [160 170]; [170 180]];

% Define broadband extraction method (see help extractBroadband_nyu)
bbmethod = 5;

% Extract broadband
for cTask = 1:length(data.task_label) 
    disp(data.task_label{cTask});
    for cRun = 1:data.nruns(cTask)
        signal = data.signal{cTask,cRun};
        broadband = extractBroadband_nyu(signal', data.fsample, bbmethod, bands);
        data.hband_sig{cTask,cRun} = broadband';
    end
end

%% PLOT: filtered time course for channel(s) of interest

task_to_plot = 'bairtemporalpattern';
runidx = 1;
taskidx = find(strcmp(task_to_plot,data.task_label));

% According to sub-beilen_acq-clinicalprojectedregions_electrodes.tsv: 
% Electrodes with lateral-occipital coverage: 
eltomatch = {'OT15'}; %, 'OT07'
% Electrodes with inferior-parietal coverage: 
%eltomatch = {'BT16', 'OT08', 'OT16'};

% Smooth & plot
figure; hold on
for ii = 1:length(eltomatch)
    % Find matching electrode numbers
    el = find(strcmp(eltomatch{ii}, data.label{taskidx,runidx}.name));
    plot(data.time{taskidx,runidx},smooth(data.hband_sig{taskidx,runidx}(el,:),128), 'LineWidth', 2);
end

% Plot event onsets on top
onsets = data.events{taskidx,runidx}.onset;
plot(onsets, zeros(length(onsets),1),'k.','MarkerSize', 25, 'LineStyle','none');
legend(eltomatch);
xlabel('time (s)');
ylabel('broadband power (60-120 Hz)');
title([data.task_label{taskidx} ' ' data.run_label{runidx}]); 

%% Epoch data (all channels, all tasks)

% Set epoch length (in seconds)
prestim = 0.5; 
poststim = 1.5;

trials = struct();
ecog_chans = find(contains(data.label{taskidx, runidx}.type, 'ECOG'));
trials.label = data.label{1,1}.name(ecog_chans);

% Determine how many samples to go back and forth to extract epoch
onset_pre = round(prestim/(1/data.fsample));
onset_post = round(poststim/(1/data.fsample));

% Extract epochs
for cTask = 1:length(data.task_label) 
    disp(data.task_label{cTask});
    for cRun = 1:data.nruns(cTask)
        [~,onsetsInx] = intersect(data.time{cTask,cRun},data.events{cTask,cRun}.onset);
        for ii = 1:length(trials.label)
            for jj = 1:length(onsetsInx)  
                
                % Broadband
                trials.broadband.(data.task_label{cTask})(ii,:,jj,cRun) = smooth(data.hband_sig{cTask,cRun}(ecog_chans(ii),onsetsInx(jj)-onset_pre:onsetsInx(jj)+onset_post),8);
                
                % Evoked
                trials.evoked.(data.task_label{cTask})(ii,:,jj,cRun) = smooth(data.signal{cTask,cRun}(ecog_chans(ii),onsetsInx(jj)-onset_pre:onsetsInx(jj)+onset_post),8);
                
                % Get event indices
                trials.events.(data.task_label{cTask})(:,cRun) = data.events{cTask,cRun}.trial_type;
            end
        end
    end
end

disp('done');

trials.time = -prestim:(1/data.fsample):poststim;
trials.fsample = data.fsample;

%% Save preprocessed data for further analysis

dataToSave = fullfile(anaDir, projectName, ['sub-' sub_label], ...
    ['sub-' sub_label '_' 'ses-' ses_label '_preproc']);

save(dataToSave, 'data', '-v7.3');

dataToSave = fullfile(anaDir, projectName, ['sub-' sub_label], ...
    ['sub-' sub_label '_' 'ses-' ses_label '_epoched']);

save(dataToSave, 'trials', '-v7.3');

