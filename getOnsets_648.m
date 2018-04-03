patientID = 648;
pth = '/Volumes/server/Projects/BAIR/ECoG/';
datapth = [pth num2str(patientID)];
stimpth = [datapth '/stimdata'];

addpath(genpath([pth '/code/']));
cd(pth);

%% specify trigger channel and task indices
trigger_channel = 138; % DC10
taskind_hrf = [1 14];
soc_taskid  = [2:13];

%% Load stimulus files
% should be 2 (hrf) + 12 (soc) stimulus files
cd(stimpth);
stimfiles = dir([num2str(patientID) '*2017*.mat']);
%stimdata = struct([]);
for ii = 1:length(stimfiles)
   stimdata(ii) = load(stimfiles(ii).name) ;    
   disp(stimdata(ii).params.loadMatrix)
   figure(ii); clf
   plot(stimdata(ii).stimulus.seqtiming, stimdata(ii).stimulus.seq, 'o-')
end
%% Load data files
cd(datapth);
% all tasks are in single data file
datafiles = dir('*.edf');
[hdr, data] = edfread(datafiles(1).name);
%[hdr.soc, data.soc] = edfread(datafiles(2).name);

%% Get Onsets

% get trigger time points from data file
sample_rate = hdr.samples(trigger_channel);
triggers = data(trigger_channel,:);
triggers = triggers / max(triggers);
[~, trigger_onsets] = findpeaks(triggers, sample_rate , 'MinPeakHeight', 0.8, 'MinPeakDistance', .5);
t = (1:size(data,2))/sample_rate; % time in seconds

% get stimulus onsets from stimulus files (concatenate all runs)

taskstimuli = []; taskindex = [];
for ii = 1:length(stimdata)
    
    stim0 = length(taskstimuli)+1;
    
    if max(taskind_hrf == ii) > 0 % hrf task
        % For hrf, we just want the time of each event. First, check these from the
        % stimulus file.
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

category_names = stimdata(2).stimulus.categories;

onsets.hrf = trigger_onsets(taskindex == 1);
onsets.soc = trigger_onsets(taskindex == 2);
stimuli.hrf = taskstimuli(taskindex == 1);
stimuli.soc = taskstimuli(taskindex == 2);

%% save 
save([patientID '_StimInfo'], 'onsets', 'stimuli', 'category_names');

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





