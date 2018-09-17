tbUse('ECoG_utils');

% SCRIPT DESCRIPTION %
% Takes BAIR data from NYU School of Medicine, gets onsets, writes out
% separate runs for each tasks, including tsv event files

%% Define paths and BIDS specs %%

patientID   = 648;
RawDataDir  = '/Volumes/server/Projects/BAIR/Data/Raw/ECoG/';
BIDSDataDir = '/Volumes/server/Projects/BAIR/Data/BIDS/';

% BIDS specs
projectName = 'visual';
sub_label   = 'som648'; 
ses_label   = 'nyuecog01';
ses_labelt1 = 'som3t01';
task_label  = {'hrfpattern',...
              'soc', ...
              'soc', ...
              'soc', ...
              'soc', ...
              'soc', ...
              'soc', ...            
              'soc', ...
              'soc', ... 
              'soc', ...              
              'soc', ...              
              'soc', ...              
              'soc', ...
              'hrfpattern'};              
run_label = {'01','01','02', '03','04', '05','06', '07','08', '09', '10','11','12', '02'};

dataWriteDir = fullfile(BIDSDataDir, projectName, sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label), 'ieeg');
stimWriteDir = fullfile(BIDSDataDir, projectName, 'stimuli');
T1WriteDir   = fullfile(BIDSDataDir, projectName, sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_labelt1), 'anat');

% define run pre and post baseline periods

prescan  = 3; % seconds
postscan = 3; % seconds

%% Read in ECoG data

% Read ECoG data
dataFiles = dir(fullfile(RawDataDir,num2str(patientID), '*.edf'));
if length(dataFiles) > 1, disp('warning: multiple datafiles found: using first one'); end

fileName = [dataFiles(1).folder filesep dataFiles(1).name];    
data = ft_read_data(fileName);
hdr = ft_read_header(fileName);
% To read in EDF data with channels with different sampling rates: data = edf2fieldtrip(fileName);

% Identify trigger channel and task indices
% for cChan = size(data,1):-1:1; figure;plot(data(cChan,:)); title([num2str(cChan) ': ' hdr.label{cChan}]);waitforbuttonpress;close; end

triggerChannel = 138; % DC10 

%% Rename and save T1 

% Locate T1
T1File = dir(fullfile(RawDataDir,num2str(patientID), 'T1.nii.gz'));
if length(T1File) > 1, disp('warning: multiple T1s found: using first one'); end
if length(T1File) < 1
    disp('warning: no T1 found!') 
else
    T1_fname = fullfile(T1WriteDir, sprintf('sub-%s_ses-%s_T1w.nii.gz', sub_label, ses_labelt1));
    [SUCCESS,MESSAGE,MESSAGEID] = copyfile(fullfile(RawDataDir,num2str(patientID), 'T1.nii.gz'), T1_fname);
end

% Note on pial file: surface recons should go in derivatives folder (see
% BIDS spec doc). If we want to include this in our data do we use our own
% freesurfer recon (which we use for electrode plotting) or the one
% provided by SOM (but why is that a .mat file)?

%% Read in electrode information

% What do we need?

% ONE coordsystem.json file 
% ONE electrodes.tsv file 

% ieeg.json file with each run (add below)
% channels.tsv file with each run (add below)

% locate and read the electrode file provided by SOM
% we want the one that matches the T1
ElecFile = dir(fullfile(RawDataDir,num2str(patientID),'*coor_T1*.txt'));
if length(ElecFile) > 1, disp('warning: multiple elec files found: using first one'); end
if length(ElecFile) < 1
    disp('warning: no coordinate file found!') 
else
    T1_fname = fullfile(T1WriteDir, sprintf('sub-%s_ses-%s_T1w.nii.gz', sub_label, ses_labelt1));
    [SUCCESS,MESSAGE,MESSAGEID] = copyfile(fullfile(RawDataDir,num2str(patientID), 'T1.nii.gz'), T1_fname);
end

D = dir([patientDir num2str(specs.pID) '/*coor_T1*.txt']);
if ~isempty(D)
    elec_file = D(1).name;
    disp(['reading ' elec_file]); % if there are multiple coor_T1 files for separate hemisphere, D(1) will always be the full list
    fid = fopen([patientDir num2str(specs.pID) '/' elec_file]); E = textscan(fid, '%s%f%f%f%s'); fclose(fid);
    elec_xyz = [E{2} E{3} E{4}]; 
    elec_labels = E{1};
    elec_types = unique(E{5});
    disp(['types of electrodes found: ' elec_types{:}]);
else
    disp('no coordinate file found - exiting');
    out = [];
    return
end

%% Get trigger time points from data file

triggers = data(triggerChannel,:);
triggers = triggers / max(triggers);
[~,trigger_onsets] = findpeaks(triggers, hdr.Fs, 'MinPeakHeight', 0.8, 'MinPeakDistance', .5);
t = ((1:hdr.nSamples)/hdr.Fs); % time in seconds

%% Read in stimulus files
% Should be 2 (hrf) + 12 (soc) stimulus files

stimDir = fullfile(RawDataDir,num2str(patientID), 'stimdata');
stimFiles = dir(fullfile(stimDir, sprintf('%s*2017*.mat', num2str(patientID))));

for ii = 1:length(stimFiles)
    fileName = [stimDir filesep stimFiles(ii).name];
    stimData(ii) = load(fileName) ;    
    disp(stimData(ii).params.loadMatrix)
    figure(ii); clf
    plot(stimData(ii).stimulus.seqtiming, stimData(ii).stimulus.seq, 'o-')
end

% CHECK: Do we have all the stimfiles?
assert(isequal(length(stimData), length(run_label)))

%% BIG LOOP ACROSS RUNS 
% This loop:
%   - gets stimulus information from stimulus files
%   - matches stimuli to trigger onsets
%   - splits data up in separate runs
% It generates:
%   - events.tsv file per run
%   - ieeg.edf file per run
%   - .mat file with stimulus info per run (goes in main 'stimuli' folder)
%   - ieeg.json file per run TO DO
%   - channels.tsv file per run TO DO

num_trials_total = 0;

for ii = 1:length(stimData)
    
    stim0 = num_trials_total+1;
    
    % Generate filename
    fname = sprintf('sub-%s_ses-%s_task-%s_run-%s', ...
            sub_label, ses_label, task_label{ii}, run_label{ii});
    fprintf('Creating edf, tsv, json and channel files for %s', fname);
    
    % Generate a json defaults struct 
    [ieeg_json, json_options] = createBIDS_ieeg_json_nyuSOM();

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
        
        % task-specific input for json file
        ieeg_json.TaskName = 'bair_hrfpattern';
        ieeg_json.TaskDescription = 'Visual textures presented at irregular intervals';
        ieeg_json.Instructions = 'Detect color change (red/green) at fixation';
        
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
        
         % task-specific input for json file
        ieeg_json.TaskName = 'bair_soc';
        ieeg_json.TaskDescription = 'Visual textures, gratings, houses faces and letter presented at various durations and inter-stimulus-intervals';
        ieeg_json.Instructions = 'Detect color change (red/green) at fixation';
        
    end
    
    % Get the onsets
    num_trials_total = num_trials_total + num_trials;
    onsets = trigger_onsets(stim0:num_trials_total);

    % Clip data from run using prescan postscan intervals
    runstart = onsets(1)-prescan; % seconds
    runstop = onsets(end)+postscan; % seconds
    [~,run_start_inx] = intersect(t,runstart); % index
    [~,run_stop_inx] = intersect(t,runstop); % index

    data_thisrun = data(:,run_start_inx:run_stop_inx);
    hdr_thisrun = hdr;
    hdr_thisrun.nSamples = size(data_thisrun,2);

    % Write out new data file
    data_fname = fullfile(dataWriteDir, sprintf('%s_ieeg.edf', fname));
    ft_write_data(data_fname, data_thisrun, 'header', hdr_thisrun, 'dataformat', 'edf');

    % Collect info for tsv files
    onset      = onsets'-onsets(1)+prescan;
    stim_file  = repmat(fname, num_trials, 1);

    % Write out tsv file 
    tsv_thisrun = table(onset, duration, ISI, trial_type, trial_name, stim_file, stim_file_index);
    tsv_fname = fullfile(dataWriteDir, sprintf('%s_events.tsv', fname));
    writetable(tsv_thisrun, tsv_fname, 'FileType','text', 'Delimiter', '\t')
    
    % Write out stimulus file
    stimfile_thisrun = stimData(ii);
    stimfile_fname = fullfile(stimWriteDir, sprintf('%s.mat', fname));
    save(stimfile_fname, '-struct', 'stimfile_thisrun', '-v7.3')
    
    % Collect info for json_ieeg file
    ieeg_json.SamplingFrequency = hdr_thisrun.Fs;
    ieeg_json.RecordingDuration = hdr_thisrun.nSamples/hdr_thisrun.Fs;
    
    % Write out json_ieeg file
    jsonfile_fname = fullfile(dataWriteDir, sprintf('%s_ieeg.json', fname));    
    jsonwrite(jsonfile_fname,ieeg_json,json_options)
    
    % CHECK: Are the onsets from the stimulus file and triggers aligned?
    figure,
    stem(trigger_onsets(stim0:stim0+num_trials-1)-trigger_onsets(stim0)); 
    hold on,
    stem(t_thisrun - t0, ':diamondr')
    xlabel('Event number'); ylabel('Time (s)');
    legend('Triggers', 'Stimulus Onsets')
end

% CHECK: Do number of triggers derived from EDF match number
% of trials from stimulus files?
assert(isequal(length(trigger_onsets), num_trials_total))



