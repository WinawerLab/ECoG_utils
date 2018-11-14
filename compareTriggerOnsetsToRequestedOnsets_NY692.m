%% NY692

%% load ECOG data
fileName = '/Volumes/server/Projects/BAIR/Data/Raw/ECoG/692/NY692_Winawer.EDF';
hdr = ft_read_header(fileName);
data = ft_read_data(fileName);
triggerChannel = 153; % DC1
t = ((1:hdr.nSamples)/hdr.Fs); % time in seconds

%% Get trigger time points from data file
triggers = data(triggerChannel,:);
triggers = triggers / max(triggers);

disp('Reading triggers from data');
[~,trigger_onsets] = findpeaks(triggers, hdr.Fs, 'MinPeakHeight', 0.8, 'MinPeakDistance', 0.5);
triggerTimes = trigger_onsets(1:226); % get triggers from just the first run (prf)

%% Plot found trigger onsets on top of trigger channel
figure('Name', 'NY692 triggers'); hold on
%plot(t,data(triggerChannel,:))
plot(t(1:120000),triggers(1:120000));
plot(triggerTimes, ones(length(triggerTimes),1)*0.9,'r.','MarkerSize', 25, 'LineStyle','none');
%stem(trigger_onsets,ones(length(trigger_onsets),1)*(max(data(triggerChannel,:))),'r');
legend({'trigger data', 'trigger onsets'}); xlabel('time (s)'); ylabel('amplitude');
title('NY 692 PRF run 1');

%% load stimulus data
load('/Volumes/server/Projects/BAIR/Data/Raw/ECoG/692/stimdata/sub-692_ses-nyuecog01_task-prf_run-1.mat');
onsetIndices = stimulus.trigSeq>0;

%% compare detected triggers, requested triggers, and flips

% get requested onsets and flip onsets
requestedTimes = stimulus.seqtiming(onsetIndices);
flipTimes = response.flip(onsetIndices);

% plot onsets 
figure('Name', 'NY692 trigger comparisons');hold on
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