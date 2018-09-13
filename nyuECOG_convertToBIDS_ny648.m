tbUse('ECoG_utils');

% SCRIPT DESCRIPTION %
% Takes BAIR data from NYU School of Medicine, gets onsets, writes out
% separate runs for each tasks, including tsv event files
% TO DO adds BIDS Metadata (json and electrode tsv files)

%% Define paths and BIDS specs %%

patientID = 648;
RawDataDir = '/Volumes/server/Projects/BAIR/Data/Raw/ECoG/';
BIDSDataDir = '/Volumes/server/Projects/BAIR/Data/BIDS/';

% BIDS specs
projectName = 'visual';
sub_label = 'som648'; 
ses_label = 'nyuecog01';
task_label = {'bairhrfpattern',...
              'bairsoc', ...
              'bairsoc', ...
              'bairsoc', ...
              'bairsoc', ...
              'bairsoc', ...
              'bairsoc', ...            
              'bairsoc', ...
              'bairsoc', ... 
              'bairsoc', ...              
              'bairsoc', ...              
              'bairsoc', ...              
              'bairsoc', ...
              'bairhrfpattern'};              
run_label = {'01','01','02', '03','04', '05','06', '07','08', '09', '10','11','12', '02'};

dataWriteDir = fullfile(BIDSDataDir, projectName, sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label), 'ieeg');
stimWriteDir = fullfile(BIDSDataDir, projectName, 'stimuli');

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
[~,trigger_onsets] = findpeaks(triggers, hdr.Fs, 'MinPeakHeight', 0.8, 'MinPeakDistance', .5);
t = ((1:hdr.nSamples)/hdr.Fs); % time in seconds

%% Read in stimulus files
% Should be 2 (hrf) + 12 (soc) stimulus files

stimDir = [RawDataDir num2str(patientID) '/stimdata/'];
stimFiles = dir([stimDir num2str(patientID) '*2017*.mat']);

for ii = 1:length(stimFiles)
    fileName = [stimDir filesep stimFiles(ii).name];
    stimData(ii) = load(fileName) ;    
    disp(stimData(ii).params.loadMatrix)
    figure(ii); clf
    plot(stimData(ii).stimulus.seqtiming, stimData(ii).stimulus.seq, 'o-')
end

%% Get stimulus information from stimulus files (concatenate all runs), match to trigger onsets

prescan  = 3; % seconds
postscan = 3; % seconds

num_trials_total = 0;

for ii = 1:length(stimData)
    
    stim0 = num_trials_total+1;
    
    % generate filename
    fname = sprintf('sub-%s_ses-%s_task-%s_run-%s', ...
            sub_label, ses_label, task_label{ii}, run_label{ii});

    if contains(task_label{ii}, 'hrf') %max(taskind_hrf == ii) > 0 % hrf task
                
        % get the timestamps
        blankIdx = mode(stimData(ii).stimulus.seq);
        blanks = stimData(ii).stimulus.seq == blankIdx;
        t0 = stimData(ii).stimulus.seqtiming(find(~blanks,1));
        t_thisrun = stimData(ii).stimulus.seqtiming(~blanks);
        
        num_trials = sum(~blanks);
        
        % task-specific input for tsv file
        duration   = ones(num_trials,1)*0.25;
        ISI = zeros(num_trials,1);
        trial_type = ones(num_trials,1);
        trial_name = repmat('HRFPATTERN', num_trials, 1);
        stim_file_index = stimData(ii).stimulus.seq(~blanks)';
        
    else % soc task
        
        % get the timestamps
        t0 = stimData(ii).stimulus.onsets(1);
        t_thisrun = stimData(ii).stimulus.onsets;
        
        num_trials = length(stimData(ii).stimulus.cat);
        
        % task-specific input for tsv file
        duration   = stimData(ii).stimulus.duration(stimData(ii).stimulus.cat)';
        ISI        = stimData(ii).stimulus.ISI(stimData(ii).stimulus.cat)';
        trial_type = stimData(ii).stimulus.cat';
        trial_name = stimData(ii).stimulus.categories(stimData(ii).stimulus.cat)';
        stim_file_index = nan(num_trials,1);
        for jj = 1:num_trials
            trial_seq = stimData(ii).stimulus.trial(jj).seq;
            stim_file_index(jj) = trial_seq(1);
        end      
    end
    
    % get the onsets
    num_trials_total = num_trials_total + num_trials;
    onsets = trigger_onsets(stim0:num_trials_total);

    % clip data from run using prescan postscan intervals
    runstart = onsets(1)-prescan; % seconds
    runstop = onsets(end)+postscan;
    [~,run_start_inx] = intersect(t,runstart); % index
    [~,run_stop_inx] = intersect(t,runstop);

    data_thisrun = data(:,run_start_inx:run_stop_inx);
    hdr_thisrun = hdr;
    hdr_thisrun.nSamples = size(data_thisrun,2);

    % write out new data file
    data_fname = fullfile(dataWriteDir, sprintf('%s_ieeg.edf', fname));
    ft_write_data(data_fname, data_thisrun, 'header', hdr_thisrun, 'dataformat', 'edf');

    % collect info for tsv files
    onset      = onsets'-onsets(1)+prescan;
    stim_file  = repmat(fname, num_trials, 1);

    % write out tsv file 
    tsv_thisrun = table(onset, duration, ISI, trial_type, trial_name, stim_file, stim_file_index);
    tsv_fname = fullfile(dataWriteDir, sprintf('%s_events.tsv', fname));
    writetable(tsv_thisrun, tsv_fname, 'FileType','text', 'Delimiter', '\t')
    
    % write out stimulus file
    stimfile_thisrun = stimData(ii);
    stimfile_fname = fullfile(stimWriteDir, sprintf('%s.mat', fname));
    save(stimfile_fname, '-struct', 'stimfile_thisrun', '-v7.3')
    
    % Plot the onsets from the stimulus file and the triggers to check whether
    % they are aligned
    figure,
    stem(trigger_onsets(stim0:stim0+num_trials-1)-trigger_onsets(stim0)); 
    hold on,
    stem(t_thisrun - t0, ':diamondr')
    xlabel('Event number'); ylabel('Time (s)');
    legend('Triggers', 'Stimulus Onsets')
end

% Check that number of triggers derived from EDF matches number of trials
% from stimulus files
assert(isequal(length(trigger_onsets), num_trials_total))

%data = edf2fieldtrip(fileName);

%% rename and save T1 

%% create required meta data files (cf BIDS started kit)
