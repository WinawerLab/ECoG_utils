
% SCRIPT DESCRIPTION %
% This script Takes BAIR data from NYU School of Medicine, gets onsets,
% writes out separate runs for each tasks, including tsv event files, and
% required BIDS metadata (coordsystem json and electrodes and channels tsv
% files). It is meant to be run cell-by-cell because some manual inputs are
% required for trigger channel selection and identifying noisy channels.

%% Define paths and BIDS specs %%

% Input paths specs
patientID   = 682;
RawDataDir  = '/Volumes/server/Projects/BAIR/Data/Raw/ECoG/';
%BIDSDataDir = '/Volumes/server/Projects/BAIR/Data/BIDS/';
BIDSDataDir = '/Users/winawerlab/matlab/git/bids-examples/';

% BIDS specs
%projectName = 'visual';
projectName = 'ieeg_visual_multimodal';

sub_label   = ['som' num2str(patientID)]; 
ses_label   = 'somecog01';
ses_labelt1 = 'som3t01';
task_label  = {'prf',...
               'prf', ...
               'hrfpattern', ...
               'spatialobject', ...
               'spatialobject', ...
               'temporalpattern', ...            
               'temporalpattern', ... 
               'spatialpattern', ...
               'spatialpattern', ...
               'spatialobject', ...
               'spatialobject', ...
               'temporalpattern', ...            
               'temporalpattern', ...
               'spatialpattern', ...
               'spatialpattern', ...
              };              
run_label = {'01','02','01','01','02','01','02','01','02','03','04','03','04','03','04'};
% NOTE: there is a 5th spatialobject logfile in stimdata, but that doesn't
% have any corresponding triggers in the recording, so we're ignoring it. 

% Output paths specs
dataWriteDir = fullfile(BIDSDataDir, projectName, sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label), 'ieeg');
stimWriteDir = fullfile(BIDSDataDir, projectName, 'stimuli');
T1WriteDir   = fullfile(BIDSDataDir, projectName, sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_labelt1), 'anat');

% Define temporal parameters
% --> These data have task onset and offset triggers.
prescan   = 0; % Segment each run with this amount before the first stimulus onset (seconds)
postscan  = 0; % Segment each run with this amount after the last stimulus onset (seconds)
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
hdr = ft_read_header(fileName);
data = ft_read_data(fileName);
% To read in EDF data with channels with different sampling rates: data = edf2fieldtrip(fileName);

%% MANUAL STEP:  Identify the trigger channel

% Define time axis (in seconds). First time point = 0 (this is assumed by
% the function we used to detect triggers below, and also in fieldtrip).
t = ((0:hdr.nSamples-1)/hdr.Fs); 

% Plot the raw voltage time course of each channel
switch makePlots 
    case 'yes'
        for cChan = size(data,1):-1:1; figure;plot(t,data(cChan,:)); title([num2str(cChan) ': ' hdr.label{cChan}]); waitforbuttonpress; close; end
end

% Write down the trigger channel (probably one labeled 'DC', see hdr.label)
triggerChannel = find(strcmp('DC1',hdr.label));

% Also write down any obviously bad channels (e.g. those with big spikes)
badChannels = [94 80]; %75:58 35 61

%% MANUAL STEP: Check channel selection

% Generate spectral plot; check command window output for outliers; 
inx_notEEGchans = find(contains(hdr.chantype, 'ecg'));
inx_DCchans = find(contains(hdr.label, 'DC'));
chansToPlot = setdiff(1:length(hdr.label),[inx_notEEGchans; inx_DCchans; badChannels']);
switch makePlots 
    case 'yes'
        [outliers] = ecog_plotChannelSpectra(data, chansToPlot,hdr);
end

% NOTE: outliers (identified as channels with mean power that is more that
% two standard deviations above or below the average across channels) should
% not be used to automatically identify bad channels, because channels with
% strong activation can have higher power on average! Instead, look at
% the time courses of those channels again to see what makes them stand out:

switch makePlots 
    case 'yes'
        for cChan = 1:length(outliers);figure; plot(t, data(outliers(cChan),:)); title([num2str(outliers(cChan)) ': ' hdr.label{outliers(cChan)}]); end
end

% If necessary, update badChannels:
morebadChannels = [75 76 1 77 78 65 68 69 70 35 58 57 49];
badChannels = [badChannels morebadChannels];

% Indicate the reason why channels were marked as bad
badChannelsDescriptions = {'clipped', 'clipped', ...
    'spikes', 'spikes', 'spikes', 'spikes', 'spikes', 'spikes', 'spikes', 'spikes', ...
    'lowfreqtransient', 'lowfreqtransient', ...
    'outlier', 'outlier', 'outlier'};

% All sEEG channels not labeled as bad will be labeled good
goodChannels = setdiff(chansToPlot,badChannels)';

% Check selection of good channels:
switch makePlots 
    case 'yes'
        % Check the powerplot of all the good channels, no leftover outliers?
        %ecog_plotChannelSpectra(data, goodChannels, hdr)

        % Check the timeseries of all the good channels, no noisy moments?
        figure;plot(t,data(goodChannels,:)); xlabel('Time (s)'); ylabel('Amplitude'); set(gca,'fontsize',16);
end

%% Get trigger time points from data file
triggers = data(triggerChannel,:);
triggers = triggers / max(triggers);

disp('Reading triggers from data');
[~,trigger_onsets] = findpeaks(triggers, hdr.Fs, 'MinPeakHeight', 0.8, 'MinPeakDistance', 0.5);

switch makePlots 
    case 'yes'
        figure;hold on
        %plot(t,data(triggerChannel,:))
        plot(t,triggers);
        plot(trigger_onsets, ones(length(trigger_onsets),1),'r.','MarkerSize', 25, 'LineStyle','none');
        %stem(trigger_onsets,ones(length(trigger_onsets),1)*(max(data(triggerChannel,:))),'r');
        legend({'trigger data', 'trigger onsets'}); xlabel('time (s)'); ylabel('amplitude');
end

% There is one trigger onset for a next run, which is (probably spatial
% object 5), for which a log file was written but no data recorded. So we
% need to remove this trigger:
trigger_onsets = trigger_onsets(1:end-1);

%% Read in stimulus files

stimDir = fullfile(RawDataDir,num2str(patientID), 'stimdata');
stimMatFiles = dir(fullfile(stimDir, sprintf('sub-*%d*.mat', patientID)));
stimTSVFiles = dir(fullfile(stimDir, sprintf('sub-*%d*.tsv', patientID)));

% CHECK: Do we have all the stimfiles?
disp('Checking whether have all the stimFiles')
assert(isequal(length(stimMatFiles), length(run_label)))
assert(isequal(length(stimTSVFiles), length(run_label)))
nRuns = length(run_label);

% NOTE: It's important to read the StimFiles in in the correct order.
% Because the logfiles are copied over, we can't assume that
% stimFiles(ii).date is correct. Instead, we sort them based on the field
% "experimentDateandTime" in params saved out by the stimulus code.

% Read the stimfiles
clear stimData;
runTimes = [];
for ii = 1:length(stimMatFiles)
    stimData(ii) = load([stimDir filesep stimMatFiles(ii).name]);
    disp(['Reading ' stimData(ii).fname])
    runTimes{ii} = stimData(ii).params.experimentDateandTime;
    switch makePlots 
        case 'yes'
            figure(ii); clf
            plot(stimData(ii).stimulus.seqtiming, stimData(ii).stimulus.trigSeq, 'o-')
    end
end

% Sort the stimfiles, count how many triggers were requested
[sortedRunTimes, runIndex] = sort(runTimes);
disp('Sorting runs in recorded order. New run order is: ')
stimData_sorted = stimData(runIndex);
requestedTriggerCount = 0;
for ii = 1:length(stimMatFiles)
    disp([stimData_sorted(ii).fname])
    num_triggers = length(find(stimData_sorted(ii).stimulus.trigSeq));
    requestedTriggerCount = requestedTriggerCount+num_triggers;
end

% CHECK: Does the number of requested triggers match the number of triggers
% that were detected in the trigger channel?
disp('Checking whether all triggers are in the data')
assert(isequal(requestedTriggerCount, length(trigger_onsets)))
stimData = stimData_sorted;

%read in tsv files later to replace onset with segmented trigger onsets
%readtable([stimDir filesep stimTSVFiles(ii).name], 'FileType');

%% Read in electrode information

% Locate and read the electrode file provided by SOM
% We want the one that matches the T1
ElecFile = dir(fullfile(RawDataDir,num2str(patientID),'*coor_T1*.txt'));
if length(ElecFile) > 1, disp('Warning: multiple elec files found: using first one'); end
if length(ElecFile) < 1
    disp('Warning: no coordinate file found!') 
else
    disp(['Reading electrode information from: ' ElecFile.name]); % if there are multiple coor_T1 files for separate hemispheres, D(1) will always be the full list
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
[elecInx, chanNames, chanTypes, chanUnits] = bidsconvert_getChannelSpecs(hdr, elec_labels, elec_types);

% Generate a channel table default (written out for each run in big loop below)
[channel_table] = createBIDS_ieeg_channels_tsv_nyuSOM(length(chanNames));
channel_table.name  = chanNames;
channel_table.type  = chanTypes;
channel_table.units = chanUnits;
channel_table.sampling_frequency = repmat(hdr.Fs,size(channel_table,1),1);

% Indicate which channel is the trigger channel in channels.tsv
channel_table.type{triggerChannel} = 'trig';

% Indicate channel statuses
%channel_table.status = repmat({'n/a'},height(channel_table),1);
%channel_table.status_description = repmat({'bad'},height(channel_table),1);
%channel_table = table2cell(channel_table);
ecogChannels = find(contains(channel_table.type,{'seeg', 'ecog'}));
for ii = 1:length(ecogChannels)
    channel_table.status{ecogChannels(ii)} = 'good';
end
for ii = 1:length(badChannels)
    channel_table.status{badChannels(ii)} = 'bad';
    if exist('badChannelsDescriptions', 'var')
        channel_table.status_description{badChannels(ii)} = badChannelsDescriptions{ii};
    end
end

%% START OF WRITING FILES %%

%% Create DATASET-SPECIFIC files %%%%%%%%%%%%%%%%%%

%   - T1w.nii.gz
%   - coordsystem.json
%   - electrodes.tsv

% %% Create T1w file

% CURRENTLY OFF BECAUSE MANUALLY EMPTIED FILES FOR bids_example

% 
% % Locate T1 provided by SoM
% T1File = dir(fullfile(RawDataDir,num2str(patientID), 'T1.nii.gz'));
% if length(T1File) > 1, disp('Warning: multiple T1s found: using first one'); end
% 
% % Copy file: rename and save in appropriate folder
% if length(T1File) < 1
%     disp('Warning: no T1 found!') 
% else
%     T1_name = fullfile(T1WriteDir, sprintf('sub-%s_ses-%s_T1w.nii.gz', sub_label, ses_labelt1));
%     disp(['Writing ' T1_name '...']);
%     copyfile(fullfile(RawDataDir,num2str(patientID), 'T1.nii.gz'), T1_name);
% end

% Note on pial file: surface recons should go in derivatives folder (see
% BIDS spec doc). If we want to include this in our data do we use our own
% freesurfer recon (which we use for electrode plotting) or the one
% provided by SOM (but why is that a .mat file)? If so, we should probably
% also provide an electrode locations file based on matched nodes
% (see run_E2N and electrode_to_nearest_node scripts in ECoG_utils)

%% Create coordsystem.json file

% Generate output name
coord_json_name = fullfile(dataWriteDir, sprintf('sub-%s_ses-%s_coordsystem.json', sub_label, ses_label));

% Get default values
[coordsystem_json, json_options] = createBIDS_ieeg_coordsystem_json_nyuSOM();

% Update values with specs for this subject
coordsystem_json.IntendedFor = fullfile(sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_labelt1), 'anat', ...
    sprintf('sub-%s_ses-%s_T1w.nii.gz', sub_label, ses_labelt1)); % this path must be specified relative to the project folder

% Write coordsystem.json file
disp(['Writing ' coord_json_name '...']);
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

% Generate output name
electrodes_tsv_name = fullfile(dataWriteDir, sprintf('sub-%s_ses-%s_electrodes.tsv', sub_label, ses_label));

% Write tsv file
disp(['Writing ' electrodes_tsv_name '...']);
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
[ieeg_json] = bidsconvert_addChannelCounts(ieeg_json, channel_table);

% If patient was collected at the OLD SOM location, overwrite the
% Manufacturers Model Name:
switch patientID
    case {645, 648, 661, 668, 674}
        ieeg_json.ManufacturersModelName = 'Natus NicoletOne C64';
end  

% Initialize trigger count
num_triggers_total = 0;

for ii = 1:nRuns
    
    stim0 = num_triggers_total+1;
    
    % Generate filename
    fname = sprintf('sub-%s_ses-%s_task-%s_run-%s', ...
            sub_label, ses_label, task_label{ii}, run_label{ii});
    fprintf('Writing eeg, events, stimuli, json and channel files for %s \n', fname);
    
    % Collect task-specific info for json file   
    if contains(task_label{ii}, 'hrf')                 
        ieeg_json.TaskName = 'hrf';
        ieeg_json.TaskDescription = 'Visual textures presented at irregular intervals';        
    elseif contains(task_label{ii}, 'prf')         
        ieeg_json.TaskName = 'prf';
        ieeg_json.TaskDescription = 'Visual bar apertures of textures sweeping across screen';	
    elseif contains(task_label{ii}, 'spatialpattern')        
        ieeg_json.TaskName = 'spatialpattern';
        ieeg_json.TaskDescription = 'Visual textures and gratings presented for brief, fixed durations';    
    elseif contains(task_label{ii}, 'temporalpattern') 
        ieeg_json.TaskName = 'temporalpattern';
        ieeg_json.TaskDescription = 'Visual textures presented for variable durations and inter-stimulus intervals';      
	elseif contains(task_label{ii}, 'spatialobject')       
        ieeg_json.TaskName = 'spatialobject';
        ieeg_json.TaskDescription = 'Houses faces and letter images presented brief, fixed durations';         
    end   
    ieeg_json.Instructions = 'Detect color change (red/green) at fixation';
    
    % Get the onsets in the data recordings for this run
    num_triggers = length(find(stimData_sorted(ii).stimulus.trigSeq));
    num_triggers_total = num_triggers_total + num_triggers;
    onsets = trigger_onsets(stim0:num_triggers_total); 
    [~,onset_indices] = intersect(t, onsets);

    % First trigger is task onset trigger:
    run_start_inx = onset_indices(1);
    % Last trigger is task offset trigger:
    run_stop_inx = onset_indices(end);
    % Stimulus onsets are those in between
    onset_indices = onset_indices(2:end-1);

    % Clip data from run using task onset and offset triggers
    data_thisrun = data(:,run_start_inx:run_stop_inx-1); % subtract one sample to make length of 'between run' data exactly 6 seconds
    hdr_thisrun = hdr;
    hdr_thisrun.nSamples = size(data_thisrun,2);
    
    % % Determine trial onset based on the triggers:
    % event_sample = (onset_indices - run_start_inx);
    
    % Determine trial onset based on the flips:   
    % Get the fliptimes for the requested triggers
	flips        = stimData(ii).response.flip(stimData(ii).stimulus.trigSeq>0); % in seconds
    % Align to first time-point. This is assumed to be same as the task
    % onset trigger on the basis of which the run is segmented.
    flips        = flips-flips(1);
    % Drop the task onset and offsets
    flips        = flips(2:end-1);
    % Convert to samples
    flip_indices = round(flips*hdr.Fs); % need to round because flip times do not always align with sample rate
    event_sample = flip_indices';
    
    % Collect info for tsv file
    events_table = stimData(ii).stimulus.tsv;
    if ~isfield(events_table, 'ISI')
        events_table.ISI = zeros(height(events_table),1);
    end
    
    % Overwrite onset with onsets of triggers
    events_table.sample       = event_sample;
    events_table.onset        = strtrim(cellstr(num2str(events_table.sample/hdr.Fs,['%.' num2str(nDecimals) 'f'])));
        
    % Update a number of other fields in events table
    events_table.stim_file    = repmat([fname '.mat'], height(events_table), 1);
    events_table.duration     = strtrim(cellstr(num2str(events_table.duration,['%.' num2str(nDecimals) 'f']))); 
    events_table.ISI          = strtrim(cellstr(num2str(events_table.ISI,['%.' num2str(nDecimals) 'f'])));    
    if contains(task_label{ii}, 'prf')
        for jj = 1:length(events_table.trial_name)
            events_table.trial_name{jj} = ['PRF-' events_table.trial_name{jj}];
        end
    end
    
    % Collect info for json_ieeg file
    ieeg_json.SamplingFrequency = hdr_thisrun.Fs;
    ieeg_json.RecordingDuration = round(hdr_thisrun.nSamples/hdr_thisrun.Fs,nDecimals);
  
%     % Write out new data file: 
%     data_thisrun =[]; % WRITE EMPTY FILES
%     data_fname = fullfile(dataWriteDir, sprintf('%s_ieeg', fname));
%     ft_write_data(data_fname, data_thisrun, 'header', hdr_thisrun, 'dataformat', 'brainvision_eeg');
%    
    % Write out events.tsv file: 
    events_fname = fullfile(dataWriteDir, sprintf('%s_events.tsv', fname));
    writetable(events_table, events_fname, 'FileType','text', 'Delimiter', '\t')
    
    %CURRENTLY OFF BECAUSE MANUALLY EMPTIED FILES FOR bids_example
%     % Write out stimulus file:
%     stimfile_thisrun = stimData(ii);
%     stimfile_fname = fullfile(stimWriteDir, sprintf('%s.mat', fname));
%     save(stimfile_fname, '-struct', 'stimfile_thisrun', '-v7.3')

    %CURRENTLY OFF BECAUSE MANUALLY EDITED HardwareFilters field (removed
    %slashes and double quotes added by jsonwrite function)
%     % Write out json_ieeg file:
%     jsonfile_fname = fullfile(dataWriteDir, sprintf('%s_ieeg.json', fname));    
%     jsonwrite(jsonfile_fname,ieeg_json,json_options)
%     
%     % Write out channels.tsv file:
%     channels_fname = fullfile(dataWriteDir, sprintf('%s_channels.tsv', fname));    
%     writetable(channel_table,channels_fname,'FileType','text','Delimiter','\t');
%     
end

% CHECK: Do number of triggers derived from EDF data file match the number
% of trials from the stimulus files?
assert(isequal(length(trigger_onsets), num_triggers_total))
disp('done');


