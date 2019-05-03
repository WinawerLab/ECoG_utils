function bidsconvert_writerunfiles(dataWriteDir, stimWriteDir, sub_label, ses_label, task_label, run_label, ...
    data, hdr, stimData, channel_table, trigger_onsets)

%% Create RUN-SPECIFIC files %%%%%%%%%%%%%%%%%%

% It generates:
%   - events.tsv
%   - ieeg.edf 
%   - .mat with stimulus info per run (goes in main 'stimuli' folder)
%   - ieeg.json 
%   - channels.tsv


%% Define temporal parameters

% USE THIS IF the data have task onset and offset triggers.
prescan   = 0; % Segment each run with this amount before the first stimulus onset (seconds)
postscan  = 0; % Segment each run with this amount after the last stimulus onset (seconds)
nDecimals = 4; % Specify temporal precision of time stamps in events files

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
if contains(sub_label, {'645', '648', '661', '668', '674'})
    ieeg_json.ManufacturersModelName = 'Natus NicoletOne C64';
end  

% Initialize trigger count
num_triggers_total = 0;
nRuns = length(run_label);

for ii = 1:nRuns
    
    stim0 = num_triggers_total+1;
    
    % Generate filename
    fname = sprintf('sub-%s_ses-%s_task-%s_run-%s', ...
            sub_label, ses_label, task_label{ii}, run_label{ii});
    fprintf('Writing eeg, events, stimuli, json and channel files for %s \n', fname);
    
    % Collect task-specific info for json file   
    if contains(task_label{ii}, 'hrf')                 
        ieeg_json.TaskName = 'bair_hrfpattern';
        ieeg_json.TaskDescription = 'Visual textures presented at irregular intervals';        
    elseif contains(task_label{ii}, 'prf')         
        ieeg_json.TaskName = 'bair_prf';
        ieeg_json.TaskDescription = 'Visual bar apertures of textures sweeping across screen';	
    elseif contains(task_label{ii}, 'spatialpattern')        
        ieeg_json.TaskName = 'bair_spatialpattern';
        ieeg_json.TaskDescription = 'Visual textures and gratings presented for brief, fixed durations';    
    elseif contains(task_label{ii}, 'temporalpattern') 
        ieeg_json.TaskName = 'bair_temporalpattern';
        ieeg_json.TaskDescription = 'Visual textures presented for variable durations and inter-stimulus intervals';      
	elseif contains(task_label{ii}, 'spatialobject')       
        ieeg_json.TaskName = 'bair_spatialobject';
        ieeg_json.TaskDescription = 'Houses faces and letter images presented brief, fixed durations';         
    end   
    ieeg_json.Instructions = 'Detect color change (red/green) at fixation';
    
    % Get the onsets in the data recordings for this run
    t = ((0:hdr.nSamples-1)/hdr.Fs); 
    num_triggers = length(find(stimData(ii).stimulus.trigSeq));
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
    events_table.event_sample = event_sample;
    events_table.onset        = strtrim(cellstr(num2str(events_table.event_sample/hdr.Fs,['%.' num2str(nDecimals) 'f'])));
        
     % Update a number of other fields in events table
    events_table.stim_file    = repmat([fname '.mat'], height(events_table), 1);
    events_table.duration     = strtrim(cellstr(num2str(events_table.duration,['%.' num2str(nDecimals) 'f']))); 
    events_table.ISI          = strtrim(cellstr(num2str(events_table.ISI,['%.' num2str(nDecimals) 'f'])));    
%     if contains(task_label{ii}, 'prf')
%         for jj = 1:length(events_table.trial_name)
%             events_table.trial_name{jj} = ['PRF' events_table.trial_name{jj}];
%         end
%     end

    % Add a task column to the events_table 
    events_table.task_name   = repmat(task_label{ii}, height(events_table), 1);
    
    % Collect info for json_ieeg file
    ieeg_json.SamplingFrequency = hdr_thisrun.Fs;
    ieeg_json.RecordingDuration = round(hdr_thisrun.nSamples/hdr_thisrun.Fs,nDecimals);
  
    % Write out new data file:    
    data_fname = fullfile(dataWriteDir, sprintf('%s_ieeg', fname));
    ft_write_data(data_fname, data_thisrun, 'header', hdr_thisrun, 'dataformat', 'brainvision_eeg');
   
    % Write out events.tsv file: 
    events_fname = fullfile(dataWriteDir, sprintf('%s_events.tsv', fname));
    writetable(events_table, events_fname, 'FileType','text', 'Delimiter', '\t')
    
    % Write out stimulus file:
    stimfile_thisrun = stimData(ii);
    stimfile_fname = fullfile(stimWriteDir, sprintf('%s.mat', fname));
    save(stimfile_fname,'-struct', 'stimfile_thisrun', '-v7.3')
      
    % Write out json_ieeg file:
    jsonfile_fname = fullfile(dataWriteDir, sprintf('%s_ieeg.json', fname));    
    jsonwrite(jsonfile_fname,ieeg_json,json_options)
    
    % Write out channels.tsv file:
    channels_fname = fullfile(dataWriteDir, sprintf('%s_channels.tsv', fname));    
    writetable(channel_table,channels_fname,'FileType','text','Delimiter','\t');
    
end

% CHECK: Do number of triggers derived from EDF data file match the number
% of trials from the stimulus files?
assert(isequal(length(trigger_onsets), num_triggers_total))
disp('done');

end