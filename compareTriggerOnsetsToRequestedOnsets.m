
% load stimulus data
load('/Volumes/server/Projects/BAIR/Data/BIDS/visual/stimuli/sub-som648_ses-nyuecog01_task-hrfpattern_run-01.mat');
onsetIndices = [false diff(stimulus.seq)<0];

% get requested onsets and flip onsets
requestedTimes = stimulus.seqtiming(onsetIndices);
flipTimes = response.flip(onsetIndices);

% load trigger data (extracted from trigger channel)
load('/Volumes/server/Projects/BAIR/Data/Raw/ECoG/648/triggersNY648.mat')
triggerTimes =  trigger_onsets(1:32); % first run

% plot onsets 
figure;hold on
stem(requestedTimes-requestedTimes(1),'g');
stem(flipTimes-flipTimes(1),'b');
stem(triggerTimes-triggerTimes(1),'r');

legend({'requested', 'flips', 'triggers'});
xlabel('Trial Number')
ylabel('Time (s)')
title('Onset per trial')
set(gca, 'FontSize', 12)

% plot intervals
figure;hold on
stem(diff(requestedTimes),'g');
stem(diff(flipTimes),'b');
stem(diff(triggerTimes),'r');

legend({'requested', 'flips', 'triggers'});
xlabel('Trial Number')
ylabel('Time (s)')
title('Interval per trial')
set(gca, 'FontSize', 12)

% plot differences in onsets
figure;hold on
stem((flipTimes-flipTimes(1))-(requestedTimes-requestedTimes(1)),'m');
stem((triggerTimes-triggerTimes(1))-(flipTimes-flipTimes(1)),'c');
stem((triggerTimes-triggerTimes(1))-(requestedTimes-requestedTimes(1)),'k');

legend({'flip minus requested', 'trigger minus flip', 'trigger minus requested'});
xlabel('Trial Number')
ylabel('Time (s)')
title('Difference in onsets')
set(gca, 'FontSize', 12)
set(gca, 'YLim', [-0.01 0.1])

% plot differences in intervals
figure;hold on
stem(diff(flipTimes)-diff(requestedTimes),'m');
stem(diff(triggerTimes)-diff(flipTimes),'c');
stem(diff(triggerTimes)-diff(requestedTimes),'k');

legend({'flip minus requested', 'trigger minus flip', 'trigger minus requested'});
xlabel('Trial Number')
ylabel('Time (s)')
title('Difference in intervals')
set(gca, 'FontSize', 12)
set(gca, 'YLim', [-0.01 0.01])
