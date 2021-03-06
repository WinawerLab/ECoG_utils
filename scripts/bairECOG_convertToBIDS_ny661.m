% SCRIPT DESCRIPTION %
% This script Takes BAIR data from NYU School of Medicine, gets onsets,
% writes out separate runs for each tasks, including tsv event files, and
% required BIDS metadata (coordsystem json and electrodes and channels tsv
% files). It is meant to be run cell-by-cell because some manual inputs are
% required for trigger channel selection and identifying noisy channels.

% Remarks:
% - Data is written in BVA format (.eeg, .vhdr, .vmrk), because of errors
% with EDF writing using fieldtrips ft_read_data and ft_write_data
% functions. Can also zip the BVA files for Flywheel (now turned off).
% - Trial onsets are based on the flips recorded by PsychToolbox, rather
% than the triggers recorded in the data (although the first trigger is
% used to detect onset of the task run). This is because we determined that
% the triggers sent over the audio cable at NYU SOM seem less accurate. 

% Questions:
% - Move the original data into a BAIR/Data/BIDS/sourcedata/ folder?
% - Should there be a json sidecar with the T1 file?

%% Define paths and BIDS specs %%

% Input paths specs
patientID   = 661;
RawDataDir  = '/Volumes/server/Projects/BAIR/Data/Raw/ECoG/';
BIDSDataDir = '/Volumes/server/Projects/BAIR/Data/BIDS/';

% BIDS specs
projectName = 'visual';
sub_label   = ['som' num2str(patientID)]; 
ses_label   = 'nyuecog01';
ses_labelt1 = 'som3t01';
acq_label   = 'clinical';
task_label  = {'prf',...
               'prf', ...
               'soc', ...
               'soc', ...
               'soc', ...
               'soc', ...
               'soc', ...            
               'soc', ...
               'soc', ... 
              'prf', ...
              'prf', ...
              'soc', ...              
              'soc', ...
              'soc'};              
run_label = {'01','02','01','02','03','04','05','06','07','03','04','08','09','10'};

% Output paths specs
dataWriteDir = fullfile(BIDSDataDir, projectName, sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label), 'ieeg');
stimWriteDir = fullfile(BIDSDataDir, projectName, 'stimuli');
T1WriteDir   = fullfile(BIDSDataDir, projectName, sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_labelt1), 'anat');

% Define temporal parameters
prescan   = 3; % Segment each run with this amount before the first stimulus onset (seconds)
postscan  = 3; % Segment each run with this amount after the last stimulus onset (seconds)
nDecimals = 4; % Specify temporal precision of time stamps in events files

% Make plots?
makePlots = 'no';

%% Create write directories

% Check whether we have the ECoG_utils repository on the path
if ~exist('createBIDS_ieeg_json_nyuSOM.m')
    tbUse ECoG_utils;
end
 
% Check whether we have write directories
if ~exist(dataWriteDir, 'dir'); mkdir(dataWriteDir);end
if ~exist(stimWriteDir, 'dir'); mkdir(stimWriteDir);end
if ~exist(T1WriteDir, 'dir'); mkdir(T1WriteDir);end

%% READ IN relevant data files %%%%%%%%%%%%%%%%%%

%   - ECoG data file
%   - electrode coordinate file provided by SoM
%   - stimulus log files generated by stimulus code

%% Read in ECoG data

% Read ECoG data
dataFiles = dir(fullfile(RawDataDir,num2str(patientID), '*.edf'));
if length(dataFiles) > 1, disp('Warning: multiple datafiles found: using first one'); end

fileName = [dataFiles(1).folder filesep dataFiles(1).name];    
disp(['Reading ' fileName '...']);
data = ft_read_data(fileName);
hdr = ft_read_header(fileName);
% To read in EDF data with channels with different sampling rates: data = edf2fieldtrip(fileName);

%% MANUAL STEP:  Identify the trigger channel

% Define time axis (in seconds). First time point = 0 (this is assumed by
% the function we used to detect triggers below, and also in fieldtrip).
t = ((0:hdr.nSamples-1)/hdr.Fs); 

% In this data set, there is a large section of recording prior to the
% start of the actual experiment, and we don't care if those parts are
% noisy, just parts with data. So let's not plot the first part:

startTime = 1700;
endTime = 3800;
[~, startSample] = min(abs(t-startTime));
[~, endSample]  = min(abs(t-endTime));

% Plot the raw voltage time course of each channel
switch makePlots 
    case 'yes'
        for cChan = size(data,1):-1:1; figure;plot(t(startSample:endSample),data(cChan,startSample:endSample)); ...
                title([num2str(cChan) ': ' hdr.label{cChan}]); waitforbuttonpress; close; end
end
 
%% Write down the trigger channel (probably one labeled 'DC', see hdr.label)
triggerChannel = 123; % DC10 

% Also write down any obviously bad channels (e.g. those with big spikes)
badChannels = [89 100 102 105:110]; 
% 100 and 102 have partly flat spectra; 89 has line noise peak on 200 Hz
% 105-110 show strong bursts of activity IN BETWEEN tasks!? Could be
% movement artifacts? Are these ECOG channels? --> checked with Hugh: these
% are 'subgaleal' electrodes that are not included in recon; exclude.

% Other observations:
% 15 has a LF drift before ~stimulus onset first task, movement?
% 47, 57, 60 have a multiple of these LF drifts throughout recording
% 69-71 have one (1) big spike during one of the PRF runs. Reason to remove?

%% Check channel selection

% Generate spectral plot; check command window output for outliers; 
switch makePlots 
    case 'yes'
        [outliers] = ecog_plotChannelSpectra(data(:,startSample:endSample), find(contains(hdr.chantype, 'eeg')), hdr);
end

% Are the outliers in the list of badChannels?
% NOTE: outliers (identified as channels with mean power that is more that
% two standard deviations above or below the average across channels) should
% not be used to automatically identify bad channels, because channels with
% strong activation can have higher power on average! Instead, look at
% the time courses of those channels again to see what makes them stand out:

switch makePlots 
    case 'yes'
        for cChan = 1:length(outliers);figure; plot(t(startSample:endSample), data(outliers(cChan),startSample:endSample)); title([num2str(outliers(cChan)) ': ' hdr.label{outliers(cChan)}]); end
end

%% All sEEG channels not labeled as bad will be labeled good
goodChannels = setdiff(find(contains(hdr.chantype, 'eeg')),badChannels)';

% Check selection of good channels:
switch makePlots 
    case 'yes'
        % Check the powerplot of all the good channels, no leftover outliers?
        %ecog_plotChannelSpectra(data(:,startSample:endSample), goodChannels, hdr)

        % Check the timeseries of all the good channels, no noisy moments?
        figure;plot(t(startSample:endSample),data(goodChannels,startSample:endSample)); xlabel('Time (s)'); ylabel('Amplitude'); set(gca,'fontsize',16);
end

% Indicate the reason why channels were marked as bad
% NOTE: do NOT use space in descriptions, readtable can't handle it!!
badChannelsDescriptions = {'spikes', 'spikes', 'linenoise', 'subgaleal', 'subgaleal' , 'subgaleal' , 'subgaleal', 'subgaleal', 'subgaleal'};
%repmat({'spikes'}, 1,length(badChannels)); 
assert(length(badChannelsDescriptions) == length(badChannels));

%% Get trigger time points from data file
triggers = data(triggerChannel,:);
triggers = triggers / max(triggers);

disp('Reading triggers from data');
%[~,trigger_onsets] = findpeaks(triggers, hdr.Fs, 'MinPeakHeight', 0.8, 'MinPeakDistance', .5);
[~,trigger_onsets] = findpeaks(triggers, hdr.Fs, 'MinPeakHeight', 0.8, 'MinPeakDistance', 0.05);

% we are looking for 4*542 PRF + 10*36 SOC = 2528 triggers
% need to adapt the minimal distance relative to NY648 because of the PRF
% triggers being closer together

switch makePlots 
    case 'yes'
        figure;hold on
        plot(t,triggers)
        stem(trigger_onsets,ones(length(trigger_onsets),1)*0.99,'r');
        legend({'trigger channel', 'detected trigger onsets'}); xlabel('time (s)'); ylabel('amplitude');
        set(gca, 'ylim', [0 1]);
end

%% Read in stimulus files

stimDir = fullfile(RawDataDir,num2str(patientID), 'stimdata');
stimFiles = dir(fullfile(stimDir, sprintf('*%s*.mat', num2str(patientID))));
clear stimData
runTimes = [];
for ii = 1:length(stimFiles)
    fileName = [stimDir filesep stimFiles(ii).name];
    stimData(ii) = load(fileName) ;    
    disp(['Reading ' stimData(ii).params.loadMatrix])
    switch makePlots 
        case 'yes'
            figure; hold on
            plot(stimData(ii).stimulus.seqtiming, stimData(ii).stimulus.seq, 'o-')
    end
    tmp = strsplit(stimFiles(ii).name,'_');
    runTimes{ii} = tmp{2}(1:end-4);
end

% CHECK: Do we have all the stimfiles?
assert(isequal(length(stimData), length(run_label)))
nRuns = length(run_label);

%% Read in electrode information

% Locate and read the electrode file provided by SOM
% We want the one that matches the T1
ElecFile = dir(fullfile(RawDataDir,num2str(patientID),'*coor_T1*.txt'));
if length(ElecFile) > 1, disp('Warning: multiple elec files found: using first one'); end
if length(ElecFile) < 1
    disp('Warning: no coordinate file found!') 
else
    disp(['Reading ' ElecFile.name]); % if there are multiple coor_T1 files for separate hemispheres, D(1) will always be the full list
    elec_fname = fullfile(ElecFile.folder,ElecFile.name);
    fid = fopen(elec_fname); E = textscan(fid, '%s%f%f%f%s'); fclose(fid);
    elec_xyz = [E{2} E{3} E{4}]; 
    elec_labels = E{1};
    elec_types = E{5};
    disp(['Types of electrodes found: ' unique(strcat(elec_types{:}))]);
    % Replace elec_types with BIDS terminology
    elec_types = replace(elec_types, 'G', 'surface');
    elec_types = replace(elec_types, 'D', 'depth');
    elec_types = replace(elec_types, 'S', 'surface');
end

% Match elec names; sort coordinates and electrode types based on matching
[elecInx, chanNames, chanTypes, chanUnits] = bidsconvert_getchannelspecs(hdr, elec_labels, elec_types);

% Generate a channel table default (written out for each run in big loop below)
[channel_table] = createBIDS_ieeg_channels_tsv_nyuSOM(length(chanNames));
channel_table.name  = chanNames;
channel_table.type  = chanTypes;
%channel_table.units = chanUnits;
channel_table.sampling_frequency = repmat(hdr.Fs,size(channel_table,1),1);

% Indicate which channel is the trigger channel in channels.tsv
channel_table.type{triggerChannel} = 'trig';

% Indicate channel statuses
%channel_table.status = cell(height(channel_table),1);
%channel_table.status_description = cell(height(channel_table),1);
nonEcogChannels = find(~contains(channel_table.type,{'seeg', 'ecog'}));
for ii = 1:length(nonEcogChannels)
    channel_table.status{nonEcogChannels(ii)} = 'bad';
end
for ii = 1:length(badChannels)
    channel_table.status{badChannels(ii)} = 'bad';
    if exist('badChannelsDescriptions', 'var')
        channel_table.status_description{badChannels(ii)} = badChannelsDescriptions{ii};
    end
end
disp('electrode information read in succesfully');
disp(channel_table);

% Indicate which channels below to which group in channels.tsv
INXNames = {'depth', 'grid' 'strip'};
INX = [];
% DEPTH electrodes:
INX{1} = find(contains(lower(channel_table.type), 'seeg'));
% GRID electrodes:
INX{2} = find(contains(lower(channel_table.type), 'ecog') & strncmp('G', channel_table.name, 1));
% ALL OTHER SURFACE electrodes
INX{3} = find(contains(lower(channel_table.type), 'ecog') & ~strncmp('G', channel_table.name, 1));
for ii = 1:length(INX)
    chan_index = INX{ii};
    fprintf('[%s] Assigning electrodes to group: %s...\n',mfilename, INXNames{ii});

    if ~isempty(chan_index)        
        channel_table.group(chan_index,:) = INXNames(ii);
    end
end

%% Create DATASET-SPECIFIC files %%%%%%%%%%%%%%%%%%

%   - T1w.nii.gz
%   - coordsystem.json
%   - electrodes.tsv

%% Create T1w file

% Locate T1 provided by SoM
T1File = dir(fullfile(RawDataDir,num2str(patientID), 'T1.nii.gz'));
if length(T1File) > 1, disp('Warning: multiple T1s found: using first one'); end

% Copy file: rename and save in appropriate folder
if length(T1File) < 1
    disp('Warning: no T1 found!') 
else
    T1_name = sprintf('sub-%s_ses-%s_acq-%s_run-01_T1w.nii.gz', sub_label, ses_labelt1, acq_label);
    T1_fname = fullfile(T1WriteDir, T1_name);
    disp(['writing ' T1_name '...']);
    copyfile(fullfile(RawDataDir,num2str(patientID), 'T1.nii.gz'), T1_fname);
    
    % Generate scans table
	filename = string(fullfile('anat', T1_name));
    acq_time = "n/a"; % We don't know when the T1 was done
    T1_scans_table = table(filename, acq_time);
    % Generate output name (hacky but we need to go one level up from anat)
    tmp = strsplit(fileparts(T1WriteDir), '/');
    T1sessionDir = strjoin(tmp,'/');
    T1_scans_name = fullfile(T1sessionDir, sprintf('sub-%s_ses-%s_scans.tsv', sub_label, ses_labelt1));
    writetable(T1_scans_table,T1_scans_name,'FileType','text','Delimiter','\t');
end

%% Create coordsystem.json file

% Generate output name
coord_json_name = fullfile(dataWriteDir, sprintf('sub-%s_ses-%s_acq-%s_coordsystem.json', sub_label, ses_label, acq_label));

% Get default values
[coordsystem_json, json_options] = createBIDS_ieeg_coordsystem_json_nyuSOM();

% Update values with specs for this subject
coordsystem_json.IntendedFor = fullfile(sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_labelt1), 'anat', ...
        T1_name); % this path must be specified relative to the project folder

% Write coordsystem.json file
disp(['writing ' coord_json_name '...']);
jsonwrite(coord_json_name,coordsystem_json,json_options);

% Note on coordystem json file: SOM also provides electrode locations in
% MNI space, as well as a list of matched anatomical regions. Do we want to
% include this information in our BIDS dataset as well? The MNI coordinates
% could be included as another coordsystem file (and associated
% electrodes.tsv file) which can be distinguished from the T1 coordinates
% by adding a [_space-<label>] to the filenames (see BIDS iEEG section
% 3.4.1). SOM does not provide an MNI aligned T1 to go along with this file
% however (unless that is the *_elec_T1.nii.gz file that they include?).

%% Create electrodes.tsv file

% Generate table with SOM defaults
n = length(elecInx);
[electrode_table] = createBIDS_ieeg_electrodes_tsv_nyuSOM(n);

% Fill with relevant information:
electrode_table.name = elec_labels(elecInx);
electrode_table.x = elec_xyz(elecInx,1);
electrode_table.y = elec_xyz(elecInx,2);
electrode_table.z = elec_xyz(elecInx,3);
electrode_table.type = elec_types(elecInx);
electrode_table.hemisphere(electrode_table.x <0) = {'L'};
electrode_table.hemisphere(electrode_table.x >0) = {'R'};

% Update electrodes.tsv to have same groups as channels.tsv
elec_inx = ecog_matchChannels(electrode_table.name, channel_table.name);
if length(find(elec_inx>0)) ~= height(electrode_table)
    error('[%s] Could not find all electrode names in channel table!', mfilename)
end
electrode_table.group = channel_table.group(elec_inx);

% Generate output name
electrodes_tsv_name = fullfile(dataWriteDir, sprintf('sub-%s_ses-%s_acq-%s_electrodes.tsv', sub_label, ses_label, acq_label));

% Write tsv file
disp(['writing ' electrodes_tsv_name '...']);
writetable(electrode_table,electrodes_tsv_name,'FileType','text','Delimiter','\t');

%% Create RUN-SPECIFIC files %%%%%%%%%%%%%%%%%%

% It generates:
%   - events.tsv
%   - ieeg.edf 
%   - .mat with stimulus info per run (goes in main 'stimuli' folder)
%   - ieeg.json 
%   - channels.tsv

%% Big loop across runs

% This loop:
%   - extracts run-specific stimulus information from stimulus files
%   - matches stimuli to trigger onsets
%   - splits data up in separate runs
%   - writes out all the run specific files required for BIDS

% Generate an _ieeg.json file default
[ieeg_json, json_options] = createBIDS_ieeg_json_nyuSOM();

% Add a count of all the channels (based on channel_table)
[ieeg_json] = bidsconvert_addchannelcounts(ieeg_json, channel_table);

% If patient was collected at the OLD SOM location, overwrite the
% Manufacturers Model Name:
switch patientID
    case {645, 648, 661, 668, 674}
        ieeg_json.ManufacturersModelName = 'Natus NicoletOne C64';
end       

% Initialize trigger count
num_triggers_total = 0;
dataFileNames = cell(nRuns,1);

for ii = 1:nRuns
    
    stim0 = num_triggers_total+1;
    
    % Generate filename
    fname = sprintf('sub-%s_ses-%s_task-%s_acq-%s_run-%s', ...
            sub_label, ses_label, task_label{ii}, acq_label, run_label{ii});
    fprintf('Writing eeg, events, stimuli, json and channel files for %s \n', fname);
    
    % Collect task-specific info for json file   
    if contains(task_label{ii}, 'prf') % prf task
        
        blankIdx = mode(stimData(ii).stimulus.seq);
        blanks = stimData(ii).stimulus.seq == blankIdx;
        num_trials = sum(~blanks);

        % task-specific input for tsv file
        duration   = stimData(ii).stimulus.tsv.duration;
        ISI        = zeros(height(stimData(ii).stimulus.tsv),1);
        trial_type = stimData(ii).stimulus.tsv.trial_type;
        trial_name = repmat('PRF',height(stimData(ii).stimulus.tsv), 1);

%        trial_name = cell(height(stimData(ii).stimulus.tsv),1);
%         for jj = 1:height(stimData(ii).stimulus.tsv)
%             trial_name{jj} = ['PRF' num2str(stimData(ii).stimulus.tsv.trial_name(jj,:))];
%         end
        stim_file_index = stimData(ii).stimulus.tsv.stim_file_index;
        
        % task-specific input for json file
        ieeg_json.TaskName = 'bair_prf';
        ieeg_json.TaskDescription = 'Visual bar apertures of textures sweeping across screen';
        
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
    
    % Get the onsets in the data recordings
    num_triggers_total = num_triggers_total + num_trials;
    onsets = trigger_onsets(stim0:num_triggers_total);
    [~,onset_indices] = intersect(t, onsets);
    
    % PRF ONLY
    switch task_label{ii}
        case 'prf'
            % first two triggers (for stim blank code 1) were not written
            % into data, so first onset_indices is actually stim 2 in
            % stimulus.tsv. We need to add additional indices to match tsv
            stimIdx = find(~blanks);
            timefromMissingTrigger = stimData(ii).stimulus.seqtiming(stimIdx(1))-stimData(ii).stimulus.seqtiming(blanks([1 2]));
            onset_indices = [onset_indices(1)-(timefromMissingTrigger*hdr.Fs)'; onset_indices];
            trig_indices = find(stimData(ii).stimulus.trigSeq>0);
            trig_indices = [1 trig_indices(1:2:end)];
            num_trials = height(stimData(ii).stimulus.tsv);
        case 'soc'
            trig_indices = stimData(ii).stimulus.trigSeq>0;
    end
    
    % Clip data from run using prescan postscan intervals
    run_start_inx = onset_indices(1)-prescan*hdr.Fs;
    run_stop_inx = onset_indices(end)+postscan*hdr.Fs;    
    data_thisrun = data(:,run_start_inx:run_stop_inx-1); % subtract one sample to make length of 'between run' data exactly 6 seconds
    hdr_thisrun = hdr;
    hdr_thisrun.nSamples = size(data_thisrun,2);
  
    % Collect info for events.tsv file 
    
    % % Determine trial onset based on the triggers:
    % event_sample = (onset_indices - run_start_inx);
    
    % Determine trial onset based on the flips:   
    % Get the fliptimes for the requested triggers
	flips        = stimData(ii).response.flip(trig_indices); % in seconds
    % Align to the flip for the first stimulus. This is assumed to be same
    % as the first stimulus trigger on the basis of which the run is
    % segmented (minus prescan period).
    flips        = flips-(flips(1)-prescan); % in seconds
    flip_indices = round(flips*hdr.Fs); % need to round because flip times do not always exactly align with sample rate
    event_sample = flip_indices';
    
    % Get remaining columns: 
    onset        = strtrim(cellstr(num2str(event_sample/hdr.Fs,['%.' num2str(nDecimals) 'f']))); 
    stim_file    = repmat([fname '.mat'], num_trials, 1);
    duration     = strtrim(cellstr(num2str(duration,['%.' num2str(nDecimals) 'f'])));
    ISI          = strtrim(cellstr(num2str(ISI,['%.' num2str(nDecimals) 'f'])));
    
    % Write out new data file  
    % % EDF format --> has problems with the SoM data (errors)
    %data_fname = fullfile(dataWriteDir, sprintf('%s_ieeg.edf', fname));
    %ft_write_data(data_fname, data_thisrun, 'header', hdr_thisrun, 'dataformat', 'edf');
    % BVA format
    data_fname = fullfile(dataWriteDir, sprintf('%s_ieeg', fname));
    ft_write_data(data_fname, data_thisrun, 'header', hdr_thisrun, 'dataformat', 'brainvision_eeg');
    % % May need to zip files for Flywheel upload (check with Gio)
    %zip(data_fname, {sprintf('%s.eeg', data_fname), sprintf('%s.vhdr', data_fname), sprintf('%s.vmrk', data_fname)});
    %delete(sprintf('%s.eeg', data_fname), sprintf('%s.vhdr', data_fname), sprintf('%s.vmrk', data_fname));
  
    % Write out tsv file 
    events_table = table(onset, duration, ISI, trial_type, trial_name, stim_file, stim_file_index, event_sample);
    
    % Add a task column to the events_table 
    events_table.task_name   = repmat(task_label{ii}, height(events_table), 1);
    
    events_fname = fullfile(dataWriteDir, sprintf('%s_events.tsv', fname));
    writetable(events_table, events_fname, 'FileType','text', 'Delimiter', '\t')
    
    % Write out stimulus file
    stimfile_thisrun = stimData(ii);
    stimfile_fname = fullfile(stimWriteDir, sprintf('%s.mat', fname));
    save(stimfile_fname, '-struct', 'stimfile_thisrun', '-v7.3')
    
    % Collect info for json_ieeg file
    ieeg_json.SamplingFrequency = hdr_thisrun.Fs;
    ieeg_json.RecordingDuration = round(hdr_thisrun.nSamples/hdr_thisrun.Fs,nDecimals);
    
    % Write out json_ieeg file
    jsonfile_fname = fullfile(dataWriteDir, sprintf('%s_ieeg.json', fname));    
    jsonwrite(jsonfile_fname,ieeg_json,json_options)
    
    % Write out channels.tsv file
    channels_tsv_fname = fullfile(dataWriteDir, sprintf('%s_channels.tsv', fname));    
    writetable(channel_table,channels_tsv_fname,'FileType','text','Delimiter','\t');
    
    % CHECK: Are the onsets from the stimulus file and triggers aligned?
    switch makePlots 
        case 'yes'
            figure; hold on;
            stem(trigger_onsets(stim0:stim0+num_trials-1)-trigger_onsets(stim0)); 
            stem(t_thisrun - t0, ':diamondr')
            xlabel('Event number'); ylabel('Time (s)');
            legend('Triggers', 'Stimulus Onsets')
    end
    
    % Collect filenames to be used for scants.tsv in
    % bidsconvert_writesessionfiles.m
    fname = fullfile('ieeg',sprintf('%s_ieeg', fname));
	dataFileNames{ii} = fname;

end

% CHECK: Do number of triggers derived from EDF data file match the number
% of trials from the stimulus files?
assert(isequal(length(trigger_onsets), num_triggers_total))
disp('done');

%% Create scans.tsv files

 % ieeg

 % Generate scans table
 acq_time = cell(length(runTimes),1);
 for ii = 1:length(runTimes)
     tt = runTimes{ii};
     tt = datestr(datenum(tt, 'yyyymmddTHHMMSS'), 'yyyy-mm-ddTHH:MM:SS');
     % shift date
     tt(1:4) = '1900';  
     acq_time{ii} = tt;
 end
 filename = dataFileNames;
 ieeg_scans_table = table(filename, acq_time);
 % Generate output name (hacky but we need to go one level up from ieeg)
 tmp = strsplit(fileparts(dataWriteDir), '/');
 ieegsessionDir = strjoin(tmp,'/');
 ieeg_scans_name = fullfile(ieegsessionDir, sprintf('sub-%s_ses-%s_scans.tsv', sub_label, ses_label));
 writetable(ieeg_scans_table,ieeg_scans_name,'FileType','text','Delimiter','\t');

