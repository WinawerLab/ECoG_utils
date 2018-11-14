% NY648 

% load trigger data (extracted from trigger channel)
%load('/Volumes/server/Projects/BAIR/Data/Raw/ECoG/648/triggersNY648.mat')

%% load ECOG data
fileName = '/Volumes/server/Projects/BAIR/Data/Raw/ECoG/648/NY648_Winawer.edf';
disp(['Reading ' fileName '...']);
data = ft_read_data(fileName);
hdr = ft_read_header(fileName);

%% Get trigger time points from data file
triggerChannel = 138; % DC10 
t = ((1:hdr.nSamples)/hdr.Fs); % time in seconds

triggers = data(triggerChannel,:);
triggers = triggers / max(triggers);

disp('Reading triggers from data');
[~,trigger_onsets] = findpeaks(triggers, hdr.Fs, 'MinPeakHeight', 0.8, 'MinPeakDistance', .5);
triggerTimes =  trigger_onsets(1:32); % get triggers from just the first run (hrf)

%% Plot found trigger onsets on top of trigger channel
figure('Name', 'NY648 triggers'); hold on
plot(t(1:205000),triggers(1:205000))
plot(triggerTimes, ones(length(triggerTimes),1)*1,'r.','MarkerSize', 25, 'LineStyle','none');
%stem(triggerTimes,ones(length(triggerTimes),1)*0.9,'r');
legend({'trigger channel', 'detected trigger onsets'}); xlabel('time (s)'); ylabel('amplitude');
%set(gca, 'ylim', [-0.1 1]);

title('NY648 HRF run 1');

%% load stimulus data
load('/Volumes/server/Projects/BAIR/Data/BIDS/visual/stimuli/sub-som648_ses-nyuecog01_task-hrfpattern_run-01.mat');
onsetIndices = [false diff(stimulus.seq)<0];

%% compare detected triggers, requested triggers, and flips

% get requested onsets and flip onsets
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
stem((triggerTimes-triggerTimes(1))-(requestedTimes-requestedTimes(1)),'k');

legend({'flip minus requested', 'trigger minus flip', 'trigger minus requested'});
xlabel('Trial Number')
ylabel('Time (ms)')
title('Difference in onsets')
set(gca, 'FontSize', 12)
%set(gca, 'YLim', [-0.01 0.1])

% plot differences in intervals
subplot(2,2,4);hold on
stem(diff(flipTimes)-diff(requestedTimes),'m');
stem(diff(triggerTimes)-diff(flipTimes),'c');
stem(diff(triggerTimes)-diff(requestedTimes),'k');

legend({'flip minus requested', 'trigger minus flip', 'trigger minus requested'});
xlabel('Trial Number')
ylabel('Time (ms)')
title('Difference in intervals')
set(gca, 'FontSize', 12)
%set(gca, 'YLim', [-0.01 0.01])


