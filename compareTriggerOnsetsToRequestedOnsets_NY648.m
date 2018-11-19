% NY648 

%% load ECOG data
fileName = '/Volumes/server/Projects/BAIR/Data/Raw/ECoG/648/NY648_Winawer.edf';
disp(['Reading ' fileName '...']);
data = ft_read_data(fileName);
hdr = ft_read_header(fileName);

%% which task?
task = 'hrf'; % soc or hrf
which_run = 1;

switch task
    case 'hrf'
        if which_run == 1
            which_trials = 1:32;
        else 
            which_trials = 465:496;
        end
        
        stim_file = sprintf('sub-som648_ses-nyuecog01_task-hrfpattern_run-%02d.mat', ...
            which_run);
    case 'soc'
        
        which_trials = 32+36*(which_run-1) + (1:36);
        stim_file = sprintf('sub-som648_ses-nyuecog01_task-soc_run-%02d.mat', ...
            which_run);
end


%% Get trigger time points from data file

triggerChannel = 138; % DC10 
%t = ((1:hdr.nSamples)/hdr.Fs); % time in seconds
t = ((0:hdr.nSamples-1)/hdr.Fs); % time in seconds

triggers = data(triggerChannel,:);
triggers = triggers / max(triggers);

disp('Reading triggers from data');
[~,trigger_onsets] = findpeaks(triggers, hdr.Fs, 'MinPeakHeight', 0.8, 'MinPeakDistance', .5);
triggerTimes =  trigger_onsets(which_trials); 
%triggerTimes = triggerTimes - 1/hdr.Fs; % shift by one sample take onset of audio wave instead of peak

%% Plot found trigger onsets on top of trigger channel
figure('Name', 'NY648 triggers'); hold on
%plot(t(1:205000),triggers(1:205000))
plot(t,triggers)
plot(triggerTimes, ones(length(triggerTimes),1)*1,'r.','MarkerSize', 25, 'LineStyle','none');
%stem(triggerTimes,ones(length(triggerTimes),1)*0.9,'r');
legend({'trigger channel', 'detected trigger onsets'}); xlabel('time (s)'); ylabel('amplitude');
%set(gca, 'ylim', [-0.1 1]);

%% load stimulus data
stim_path = '/Volumes/server/Projects/BAIR/Data/BIDS/visual/stimuli/';
load(fullfile(stim_path, stim_file));

onsetIndices = [false diff(stimulus.seq)<0];

if strcmp(task, 'soc')
   tmp = find(onsetIndices);
   duplicates = find(diff([false stimulus.seq(tmp)]) == 0);
   onsetIndices(tmp(duplicates))=0;
end

%% compare detected triggers, requested triggers, and flips

% get requested onsets and flip onsets
triggerTimes =  trigger_onsets(which_trials); % get triggers from just the first run (hrf)
requestedTimes = stimulus.seqtiming(onsetIndices);
flipTimes = response.flip(onsetIndices);

% plot onsets 
figure('Name', 'NY648 trigger comparisons');hold on
subplot(2,2,1);hold on
stem(requestedTimes-requestedTimes(1),'g');
stem(flipTimes-flipTimes(1),'b');
stem(triggerTimes-triggerTimes(1),'r');

legend({'requested', 'flips', 'triggers'});
xlabel('Trial Number')
ylabel('Time (s)')
title('Onset per trial')
set(gca, 'FontSize', 12)

% plot intervals
subplot(2,2,2);hold on
stem(diff(requestedTimes),'g');
stem(diff(flipTimes),'b');
stem(diff(triggerTimes),'r');

legend({'requested', 'flips', 'triggers'});
xlabel('Trial Number')
ylabel('Time (s)')
title('Interval per trial')
set(gca, 'FontSize', 12)

% convert to ms
flipTimes = flipTimes * 1000;
triggerTimes = triggerTimes * 1000;
requestedTimes = requestedTimes * 1000;

% plot differences in onsets
subplot(2,2,3);hold on
stem((flipTimes-flipTimes(1))-(requestedTimes-requestedTimes(1)),'m');
stem((triggerTimes-triggerTimes(1))-(flipTimes-flipTimes(1)),'c');
%stem((triggerTimes-triggerTimes(1))-(requestedTimes-requestedTimes(1)),'k');

legend({'flip minus requested', 'trigger minus flip'});
xlabel('Trial Number')
ylabel('Time (ms)')
title('Difference in onsets')
set(gca, 'FontSize', 12)
%set(gca, 'YLim', [-0.01 0.1])

% plot differences in intervals
subplot(2,2,4);hold on
stem(diff(flipTimes)-diff(requestedTimes),'m');
stem(diff(triggerTimes)-diff(flipTimes),'c');
%stem(diff(triggerTimes)-diff(requestedTimes),'k');

legend({'flip minus requested', 'trigger minus flip'});
xlabel('Trial Number')
ylabel('Time (ms)')
title('Difference in intervals')
set(gca, 'FontSize', 12)
%set(gca, 'YLim', [-0.01 0.01])

%% compare accumulated error in flips and triggers between first stimulus of first run and last stimulus of last run
stim_file1 = 'sub-som648_ses-nyuecog01_task-hrfpattern_run-01.mat';
stim_file2 = 'sub-som648_ses-nyuecog01_task-hrfpattern_run-02.mat';

s1 = load(fullfile(stim_path, stim_file1));
s2 = load(fullfile(stim_path, stim_file2));

onsetIndices1 = [false diff(s1.stimulus.seq)<0];
onsetIndices2 = [false diff(s2.stimulus.seq)<0];


flips1 = s1.response.flip(onsetIndices1);
flips2 = s2.response.flip(onsetIndices2);

t0 = s1.response.flip(1);
t1 = flips1(1);
tf = flips2(end);

% accumulated time according to screen flips (PTB)
accumFlips = tf-t1;

% accumulated time according to triggers
accumTrigs = trigger_onsets(end)-trigger_onsets(1);

% accumulated mismatch in ms over the whole set of experiments
disp((accumTrigs-accumFlips)*1000)

% mismatch (in ms) per sample (ie per 1/512 s)
disp((accumTrigs-accumFlips)*1000 / (accumTrigs*512))

% mismatch (in ms) per second
disp((accumTrigs-accumFlips)*1000 / (accumTrigs))