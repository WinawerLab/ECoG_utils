patientID = 661;
pth = '/Volumes/server/Projects/BAIR/ECoG/';
datapth = [pth num2str(patientID)];
stimpth = [datapth '/stimdata'];

addpath(genpath([pth '/code/']));
cd(pth);

%% specify trigger channel and task indices
trigger_channel = 123; % DC11

% The notes from the recording session supplied by Lora, the first two
% 'blocks' received no triggers; this would correspond to our first two PRF
% runs. However, the trigger channel suggests that these runs WERE
% recorded; instead not all of the spatiotemporal runs were recorded...?

% IF going with trigger channel appearance:
taskind_prf = [1 2 10 11];
taskind_soc = [3:9 12:14]; 
% PROBLEM:  Now not sure which stimfiles correspond to runs 3:9
% SOLUTION: Compared the distance in writing time of stimdata logfiles
% between runs with the distances of the last onset of the stimuli in each
% run: concluded that ny661_20180211T114254 and ny661_20180211T114312 were
% NOT written into the edf file. So these runs are EXCLUDED from analysis.

%% Load stimulus files
cd(stimpth);
stimfiles = dir(['*' num2str(patientID) '*2018*.mat']);
for ii = 1:length(stimfiles)
   stimdata(ii) = load(stimfiles(ii).name) ;    
   disp(stimdata(ii).params.loadMatrix)
   figure(ii); clf
   plot(stimdata(ii).stimulus.seqtiming, stimdata(ii).stimulus.seq, 'o-')
end

%% Load data files

% FieldTrip is slower but does give time axis
dataName = fullfile(datapth, sprintf('NY%d_BairPilot.edf', patientID));

cfg            = [];
cfg.dataset    = dataName;
cfg.continuous = 'yes';
cfg.channel    = 'all';
data           = ft_preprocessing(cfg);

%% Get Onsets

% we are looking for 4*542 PRF + 10*36 SOC = 2528 triggers
% need to adapt the minimal distance relative to NY648 because of the PRF
% triggers being closer together
sample_rate = data.fsample;
triggers = data.trial{1}(trigger_channel,:);
triggers = triggers / max(triggers);
[~, trigger_onsets] = findpeaks(triggers, sample_rate , 'MinPeakHeight', 0.8, 'MinPeakDistance', 0.05);

figure;plot(data.time{1},triggers)
hold on
stem(trigger_onsets, ones(length(trigger_onsets),1), '-or');
axis([1400 4000 0 1.1]);

disp(['found ' num2str(length(trigger_onsets)) ' triggers']);

%% Get stimulus onsets from stimulus files (concatenate all runs); compare

taskstimuli = []; taskindex = []; onsets = []; stimuli = [];

for ii = 1:length(stimdata)
    
    stim0 = length(taskstimuli)+1;    
    if max(taskind_prf == ii) > 0 % prf task
        blankIdx = mode(stimdata(ii).stimulus.seq);
        blanks = stimdata(ii).stimulus.seq == blankIdx;
        num_trials = sum(~blanks);
        taskindex = [taskindex ones(1,num_trials)];
        taskstimuli = [taskstimuli 99*ones(size(find(~blanks)))];
        t0 = stimdata(ii).stimulus.seqtiming(find(~blanks,1));
        t = stimdata(ii).stimulus.seqtiming(~blanks);
    else % soc task
        taskstimuli = [taskstimuli stimdata(ii).stimulus.cat];
        num_trials = length(stimdata(ii).stimulus.cat);
        taskindex = [taskindex ones(1,num_trials)*2];
        t0 = stimdata(ii).stimulus.onsets(1);
        t = stimdata(ii).stimulus.onsets;
    end
    
    % Plot the onsets from the stimulus file and the triggers to check whether
    % they are aligned
    figure,
    stem(trigger_onsets(stim0:stim0+num_trials-1)-trigger_onsets(stim0)); 
    hold on,
    stem(t - t0, ':diamondr')
    xlabel('Event number'); ylabel('Time (s)');
    legend('Triggers', 'Stimulus Onsets')
end

% Check that number of triggers derived from EDF matches number of trials
% from stimulus files
assert(isequal(length(trigger_onsets), length(taskstimuli)))

category_names = stimdata(3).stimulus.categories;

onsets.prf = trigger_onsets(taskindex == 1);
onsets.soc = trigger_onsets(taskindex == 2);
stimuli.prf = taskstimuli(taskindex == 1);
stimuli.soc = taskstimuli(taskindex == 2);

%% save 
save(sprintf('NY%d_StimInfo.mat', patientID), 'onsets', 'stimuli', 'category_names');

%% to determine which category was shown in hrf task

% hrf stim
figure;imshow(stimdata(1).stimulus.images(:,:,3))

% soc stim
figure;hold on;
for cIm = 1:36
    subplot(6,6,cIm);
    imshow(stimdata(2).stimulus.im_cell{cIm}(:,:,1));
end

% HRF stim == soc stim category 10 (SPARSITY_2)? Or 11 (SPARSITY_3)? 
% Not sure!





