pth = '/Volumes/server/Projects/BAIR/BIDS_example/visualFullSet/sub-som645/ses-somIEMU01/ieeg';
cd(pth);

%% Load stimulus files
% should be 4 stimulus files
stimfiles = dir('2017*.mat');
for ii = 1:length(stimfiles)
   results(ii) = load(stimfiles(ii).name) ;    
   disp(results(ii).params.loadMatrix)
   figure(ii); clf
   plot(results(ii).stimulus.seqtiming, results(ii).stimulus.seq, 'o-')
end

%% Load data files
datafiles = dir('*.edf');
[hdr.hrf, data.hrf] = edfread(datafiles(1).name);
[hdr.soc, data.soc] = edfread(datafiles(2).name);

%% HRF Onsets

% For hrf, we just want the time of each event. First, check these from the
% stimulus file.
blankIdx = mode(results(1).stimulus.seq);
blanks   = results(1).stimulus.seq == blankIdx;
num_trials = sum(~blanks);

trigger_channel = 117;
sample_rate = hdr.hrf.samples(trigger_channel);
t.hrf = (1:size(data.hrf,2))/sample_rate; % time in seconds

triggers = data.hrf(trigger_channel,:);
triggers = triggers / max(triggers);
[~, onsets.hrf] = findpeaks(triggers, sample_rate , 'MinPeakHeight', 0.8, 'MinPeakDistance', .5);
 stimuli.hrf = 99*ones(size(onsets.hrf));
% Check that number of triggers derived from EDF matches number of trials
% from stimulus files
assert(isequal(length(onsets.hrf), num_trials))

% Plot the onsets from the stimulus file and the triggers to check whether
% they are aligned
figure,
stem(onsets.hrf-onsets.hrf(1)); 
xlabel('Event number'); ylabel('Time (s)');

hold on,
t0 = results(1).stimulus.seqtiming(find(~blanks,1));
stem(results(1).stimulus.seqtiming(~blanks) - t0, ':diamondr')

legend('Triggers', 'Stimulus Onsets')

%% Spatiotemporal (SOC) onsets

category_names = results(2).stimulus.categories;
stimuli.soc = [results(2).stimulus.cat results(3).stimulus.cat results(4).stimulus.cat];
triggers = data.soc(trigger_channel,:);
triggers = triggers / max(triggers);
[~, onsets.soc] = findpeaks(triggers, sample_rate , 'MinPeakHeight', 0.8, 'MinPeakDistance', .5);
 
