%% NY692

%% load ECOG data
fileName = '/Volumes/server/Projects/BAIR/Data/Raw/ECoG/692/NY692_Winawer.EDF';
hdr = ft_read_header(fileName);
data = ft_read_data(fileName);
triggerChannel = 153; % DC1
%t = ((1:hdr.nSamples)/hdr.Fs); % time in seconds
t = ((0:hdr.nSamples-1)/hdr.Fs); % time in seconds

%% Get trigger time points from data file
triggers = data(triggerChannel,:);
triggers = triggers / max(triggers);

disp('Reading triggers from data');
[~,trigger_onsets] = findpeaks(triggers, hdr.Fs, 'MinPeakHeight', 0.8, 'MinPeakDistance', 0.5);
%[pks, locs] = findpeaks(triggers,'MinPeakHeight', 0.8, 'MinPeakDistance', 0.5);

%% which task?
task = 'prf'; % prf, spatialpattern, spatialobject, temporalpattern
which_run = 1; % 1 or 2 (ignoring runs 3 and 4 for non-prf runs)

switch task
    case 'prf'
        which_trials = 226*(which_run-1) + (1:226); %(1:226)*which_run;
        stim_file = sprintf('sub-som692_ses-nyuecog01_task-prf_run-%02d.mat', ...
            which_run);
    case 'spatialpattern'
        which_trials = 452+38*(which_run-1) + (1:38);
        stim_file = sprintf('sub-som692_ses-nyuecog01_task-spatialpattern_run-%02d.mat', ...
            which_run);
	case 'spatialobject'
        which_trials = 452+2*38+38*(which_run-1) + (1:38);
        stim_file = sprintf('sub-som692_ses-nyuecog01_task-spatialobject_run-%02d.mat', ...
            which_run);
	case 'temporalpattern'
        which_trials = 452+4*38+38*(which_run-1) + (1:38);
        stim_file = sprintf('sub-som692_ses-nyuecog01_task-temporalpattern_run-%02d.mat', ...
            which_run);
end

triggerTimes =  trigger_onsets(which_trials); 

%% Plot found trigger onsets on top of trigger channel
figure('Name', 'NY692 triggers'); hold on
%plot(t,data(triggerChannel,:))
%plot(t(1:120000),triggers(1:120000));
plot(t,triggers);
plot(triggerTimes, ones(length(triggerTimes),1)*1,'r.','MarkerSize', 25, 'LineStyle','none');
%stem(trigger_onsets,ones(length(trigger_onsets),1)*(max(data(triggerChannel,:))),'r');
legend({'trigger data', 'trigger onsets'}); xlabel('time (s)'); ylabel('amplitude');

%% load stimulus data
stim_path = '/Volumes/server/Projects/BAIR/Data/BIDS/visual/stimuli/';
load(fullfile(stim_path, stim_file));
onsetIndices = stimulus.trigSeq>0;

%% compare detected triggers, requested triggers, and flips

% get requested onsets and flip onsets
requestedTimes = stimulus.seqtiming(onsetIndices);
flipTimes = response.flip(onsetIndices);

% plot onsets 
figure('Name', 'NY692 trigger comparisons');hold on
subplot(2,2,1);hold on
%stem(requestedTimes-requestedTimes(1),'g');
stem(flipTimes-flipTimes(1),'b');
stem(triggerTimes-triggerTimes(1),'r');

%legend({'requested', 'flips', 'triggers'});
legend({'flips', 'triggers'});
xlabel('Trial Number')
ylabel('Time (s)')
title('Onset per trial')
set(gca, 'FontSize', 12)

% plot intervals
subplot(2,2,2);hold on
%stem(diff(requestedTimes),'g');
stem(diff(flipTimes),'b');
stem(diff(triggerTimes),'r');

%legend({'requested', 'flips', 'triggers'});
legend({'flips', 'triggers'});
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
%stem((flipTimes-flipTimes(1))-(requestedTimes-requestedTimes(1)),'m');
stem((triggerTimes-triggerTimes(1))-(flipTimes-flipTimes(1)),'c');
%stem((triggerTimes-triggerTimes(1))-(requestedTimes-requestedTimes(1)),'k');

%legend({'flip minus requested', 'trigger minus flip'});
legend({'trigger minus flip'});
xlabel('Trial Number')
ylabel('Time (ms)')
title('Difference in onsets')
set(gca, 'FontSize', 12)
%set(gca, 'YLim', [-0.01 0.1])

% plot differences in intervals
subplot(2,2,4);hold on
%stem(diff(flipTimes)-diff(requestedTimes),'m');
stem(diff(triggerTimes)-diff(flipTimes),'c');
%stem(diff(triggerTimes)-diff(requestedTimes),'k');

%legend({'flip minus requested', 'trigger minus flip'});
legend({'trigger minus flip'});

xlabel('Trial Number')
ylabel('Time (ms)')
title('Difference in intervals')
set(gca, 'FontSize', 12)
%set(gca, 'YLim', [-0.01 0.01])

%% compare accumulated error in flips and triggers between first stimulus of first run and last stimulus of last run

stim_path = '/Volumes/server/Projects/BAIR/Data/BIDS/visual/stimuli/';

stim_file1 = 'sub-som692_ses-nyuecog01_task-prf_run-01.mat';
stim_file2 = 'sub-som692_ses-nyuecog01_task-temporalpattern_run-04.mat';

s1 = load(fullfile(stim_path, stim_file1));
s2 = load(fullfile(stim_path, stim_file2));

%onsetIndices1 = [false diff(s1.stimulus.seq)<0];
%onsetIndices2 = [false diff(s2.stimulus.seq)<0];
onsetIndices1 = s1.stimulus.trigSeq>0;
onsetIndices2 = s2.stimulus.trigSeq>0;

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