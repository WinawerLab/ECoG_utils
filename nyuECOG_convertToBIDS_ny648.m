tbUse('ECoG_utils');

% SCRIPT DESCRIPTION %
% Takes BAIR data from NYU School of Medicine, gets onsets, splits runs
% Adds BIDS Metadata

%% Define paths and BIDS specs %%

patientID = 648;
RawDataDir = '/Volumes/server/Projects/BAIR/Data/Raw/ECoG/';
BIDSDataDir = '/Volumes/server/Projects/BAIR/Data/BIDS/';

% BIDS specs
projectName = 'visual';
sub_label = 'som648'; 
ses_label = 'nyuECOG01';
task_label = {'hrf'; 'soc'};
run_label = {[1 2]; [1:12]}; %run_label = '01';

%% Read in ECoG data

% Read ECoG data
dataFiles = dir([RawDataDir num2str(patientID) filesep '*.edf']);
if length(dataFiles) > 1, disp('warning: multiple datafiles found: using first one'); end

fileName = [dataFiles(1).folder filesep dataFiles(1).name];    
data = ft_read_data(fileName);
hdr = ft_read_header(fileName);
% to read in EDF data with channels with different sampling rates: data = edf2fieldtrip(fileName);

% Identify trigger channel and task indices
% for cChan = size(data,1):-1:1; figure;plot(data(cChan,:)); title([num2str(cChan) ': ' hdr.label{cChan}]);waitforbuttonpress;close; end

triggerChannel = 138; % DC10 

%% Get trigger time points from data file

triggers = data(triggerChannel,:);
triggers = triggers / max(triggers);
[~, trigger_onsets] = findpeaks(triggers, hdr.Fs, 'MinPeakHeight', 0.8, 'MinPeakDistance', .5);
t = ((1:hdr.nSamples)/hdr.Fs); % time in seconds

%% Read in stimulus files
% Should be 2 (hrf) + 12 (soc) stimulus files

stimDir = [RawDataDir num2str(patientID) '/stimdata/'];
stimFiles = dir([stimDir num2str(patientID) '*2017*.mat']);

%assert(length(stimFiles) == numel(run_label{:}))

for ii = 1:length(stimFiles)
   fileName = [stimDir filesep stimFiles(ii).name];
   stimData(ii) = load(fileName) ;    
   disp(stimData(ii).params.loadMatrix)
   figure(ii); clf
   plot(stimData(ii).stimulus.seqtiming, stimData(ii).stimulus.seq, 'o-')
end

%% Get stimulus onsets from stimulus files (concatenate all runs)

taskstimuli = []; taskindex = [];
for ii = 1:length(stimData)
    
    stim0 = length(taskstimuli)+1;
    
    if max(taskind_hrf == ii) > 0 % hrf task
        % For hrf, we just want the time of each event. First, check these from the
        % stimulus file.
        blankIdx = mode(stimData(ii).stimulus.seq);
        blanks = stimData(ii).stimulus.seq == blankIdx;
        num_trials = sum(~blanks);
        taskindex = [taskindex ones(1,num_trials)];
        taskstimuli = [taskstimuli 99*ones(size(find(~blanks)))];
        t0 = stimData(ii).stimulus.seqtiming(find(~blanks,1));
        t = stimData(ii).stimulus.seqtiming(~blanks);
    else % soc task
        taskstimuli = [taskstimuli stimData(ii).stimulus.cat];
        num_trials = length(stimData(ii).stimulus.cat);
        taskindex = [taskindex ones(1,num_trials)*2];
        t0 = stimData(ii).stimulus.onsets(1);
        t = stimData(ii).stimulus.onsets;
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

category_names = stimData(2).stimulus.categories;

onsets.hrf = trigger_onsets(taskindex == 1);
onsets.soc = trigger_onsets(taskindex == 2);
stimuli.hrf = taskstimuli(taskindex == 1);
stimuli.soc = taskstimuli(taskindex == 2);

%ft_write_data('test.edf', data, 'header', hdr, 'dataformat', 'edf');

%% Reformat data files into separate runs / tasks    


%taskind_hrf = [1 14];
%soc_taskid  = [2:13];

%data = edf2fieldtrip(fileName);


 
%% Load stimulus files
% should be 2 (hrf) + 12 (soc) stimulus files
cd(stimpth);
stimfiles = dir([num2str(patientID) '*2017*.mat']);
%stimdata = struct([]);
for ii = 1:length(stimfiles)
   stimData(ii) = load(stimfiles(ii).name) ;    
   disp(stimData(ii).params.loadMatrix)
   figure(ii); clf
   plot(stimData(ii).stimulus.seqtiming, stimData(ii).stimulus.seq, 'o-')
end
%% Load data files
cd(datapth);

%[hdr.soc, data.soc] = edfread(datafiles(2).name);

% reformat

% save under appropriate subject/session/task/run labels


%% rename and save T1 

%% create required meta data files (cf BIDS started kit)
