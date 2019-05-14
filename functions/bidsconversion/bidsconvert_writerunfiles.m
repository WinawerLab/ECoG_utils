function bidsconvert_writerunfiles(dataWriteDir, stimWriteDir, sub_label, ses_label, task_label, run_label, ...
    data, hdr, stimData, channel_table, trigger_onsets)

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
if contains(sub_label, {'645', '648', '661', '668', '674'})
    ieeg_json.ManufacturersModelName = 'Natus NicoletOne C64';
end  

% Initialize trigger count
num_triggers_total = 0;
nRuns = length(run_label);
segmentOnFlips = 0;% Make this an input variable?
nDecimals = 4; % Specify temporal precision of time stamps in events files

for ii = 1:nRuns
    
    stim0 = num_triggers_total+1;
    
    % Generate filename
    fname = sprintf('sub-%s_ses-%s_task-%s_run-%s', ...
            sub_label, ses_label, task_label{ii}, run_label{ii});
    fprintf('[%s] Writing eeg, events, stimuli, json and channel files for %s \n', mfilename, fname);
    
   
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
    
    % Fix some problems in older stimulus files 
    if contains(task_label{ii}, 'hrf')   
        if length(find(stimData(ii).stimulus.trigSeq)) ~=32 && max(stimData(ii).stimulus.trigSeq) ~= 256
            % There are additional triggers that we should not use
            fprintf('[%s] Warning: hrf stimulus file has incorrect number of triggers \n', mfilename)
            fprintf('[%s] Removing directly adjacent triggers from stimfile \n', mfilename)
            requestedTriggerInx = find(stimData(ii).stimulus.trigSeq);
            requestedTriggerInx_diff = [nan diff(requestedTriggerInx)];
            stimData(ii).stimulus.trigSeq(requestedTriggerInx(requestedTriggerInx_diff == 1)) = 0;
        end
    end
    if contains(task_label{ii}, 'prf')
        if length(find(stimData(ii).stimulus.trigSeq)) ~=224 && max(stimData(ii).stimulus.trigSeq) ~= 256
            % Old prf stim files had too many triggers + missing triggers
            fprintf('[%s] Warning: prf stimulus file has incorrect number of triggers \n', mfilename)
            [newTrigSeq] = bidsconvert_fixprftrigseq(stimData(ii).stimulus);
            stimData(ii).stimulus.trigSeq = newTrigSeq;
            % Overwrite trial names to prevent concatenation problems in
            % preprocessing
            stimData(ii).stimulus.tsv.trial_name = [repmat('PRF', [height(stimData(ii).stimulus.tsv),1]) num2str(stimData(ii).stimulus.tsv.trial_name)];
        end
    end
    
    % Get the onsets in the data recordings for this run
    t = ((0:hdr.nSamples-1)/hdr.Fs); 
    num_triggers = length(find(stimData(ii).stimulus.trigSeq));
    num_triggers_total = num_triggers_total + num_triggers;
    onsets = trigger_onsets(stim0:num_triggers_total); 
    [~,onset_indices] = intersect(t, onsets);
   
    % Check if there are onset/offset triggers
    if max(stimData(ii).stimulus.trigSeq) == 256
        % THESE DATA HAVE ONSET/OFFSET TRIGGERS
        hasOnOffTriggers = 1;
        % First trigger is task onset trigger:
        run_start_inx = onset_indices(1);
        % Last trigger is task offset trigger:
        run_stop_inx = onset_indices(end);
        % Stimulus onsets are those in between
        onset_indices = onset_indices(2:end-1);
    else
        % THESE DATA DO NOT HAVE ONSET/OFFSET TRIGGERS
        hasOnOffTriggers = 0;
        % Segment 3 seconds before first onset
        run_start_inx = onset_indices(1)-3*hdr.Fs;
        % Segment 3 seconds after last onset
        run_stop_inx = onset_indices(end)+3*hdr.Fs;             
    end
    
    % Clip data from run using task onset and offset triggers
    data_thisrun = data(:,run_start_inx:run_stop_inx-1); % subtract one sample to make length of 'between run' data exactly 6 seconds
    hdr_thisrun = hdr;
    hdr_thisrun.nSamples = size(data_thisrun,2);
    
    if segmentOnFlips    
        % Determine trial onset based on the flips:   
        % Get the fliptimes for the requested triggers
        flips        = stimData(ii).response.flip(stimData(ii).stimulus.trigSeq>0); % in seconds
        % Align to first flip
        flips        = flips-flips(1);
        if hasOnOffTriggers
            % First flip is the same as the onset trigger on the basis of
            % which the run is segmented. Drop the task onset and offsets
            flips = flips(2:end-1);
        else
            % First flip is aligned with the first stimulus onset trigger,
            % which we position 3 seconds after start of the run:
            flips = flips+3;
        end
        % Convert to samples
        flip_indices = round(flips*hdr.Fs); % need to round because flip times do not always align with sample rate
        event_sample = flip_indices';
    else
        % Determine trial onset based on the triggers:
        event_sample = (onset_indices - run_start_inx);
    end
    
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
fprintf('[%s] BIDS conversion done\n', mfilename)

end