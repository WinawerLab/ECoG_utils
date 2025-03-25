function [dataFileNames] = bidsconvert_writerunfiles(dataWriteDir, stimWriteDir, ...
    sub_label, ses_label, task_label, acq_label, run_label, ...
    data, hdr, stimData, channel_table, trigger_onsets, segmentOnFlips)

%% Create RUN-SPECIFIC files %%%%%%%%%%%%%%%%%%

% It generates:
%   - events.tsv
%   - ieeg.edf 
%   - .mat with stimulus info per run (goes in main 'stimuli' folder)
%   - ieeg.json 
%   - channels.tsv
%
% Should be run BEFORE generating the session specific files (so we can generate
% the scans.tsv files based on the existing dataFileNames)

if ~exist('segmentOnFlips', 'var') || isempty(segmentOnFlips)
    segmentOnFlips = 1;
end

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

% Initialize trigger count
num_triggers_total = 0;
nRuns = length(run_label);
nDecimals = 4; % Specify temporal precision of time stamps in events files
dataFileNames = cell(nRuns,1);

for ii = 1:nRuns
    
    stim0 = num_triggers_total+1;
    
    % Check the sensory domain; assume visual (for older datasets)
    if isfield(stimData(ii).params, 'sensoryDomain')
        sensoryDomain = lower(stimData(ii).params.sensoryDomain);
    else
        sensoryDomain = 'visual';
    end
       
    % Collect task-specific info for json file   
    switch sensoryDomain
        
        case 'visual'
            % VISUAL TASKS
            if contains(task_label{ii}, 'hrf')                 
                ieeg_json.TaskName = 'hrfpattern';
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
                ieeg_json.TaskDescription = 'Houses, faces and letter images presented for brief, fixed durations';
            end   
            ieeg_json.Instructions = 'Detect color change (red/green) at fixation';
            
        case 'motor'
            % MOTOR TASKS
            if contains(task_label{ii}, 'boldhand')                 
                ieeg_json.TaskName = 'boldhand';
                ieeg_json.TaskDescription = 'Hand clenching upon upon a visual cue (fixation color change) presented at irregular intervals';  
                ieeg_json.Instructions = 'When fixation dot color changes, clench hand into a fist';
            elseif contains(task_label{ii}, 'boldsat') 
                ieeg_json.TaskName = 'boldsat';
                ieeg_json.TaskDescription = 'Hand clenching upon a visual cue (fixation color change) presented a regular intervals at one of four frequencies (boldsat 1-4)';      
                ieeg_json.Instructions = 'When fixation dot color changes, clench hand into a fist';
            elseif contains(task_label{ii}, 'fingermappingleft')         
                ieeg_json.TaskName = 'fingermapping_left';
                ieeg_json.TaskDescription = 'Moving one finger at a time based on a visual cue indicating which finger to move (left hand)';	
                ieeg_json.Instructions = 'When one of the 5 cues on the screen changes from white to black, bend the associated finger';
                task_label{ii} = 'fingermapping';
            elseif contains(task_label{ii}, 'fingermappingright')         
                ieeg_json.TaskName = 'fingermapping_right';
                ieeg_json.TaskDescription = 'Moving one finger at a time based on a visual cue indicating which finger to move (right hand)';	
                ieeg_json.Instructions = 'When one of the 5 cues on the screen changes from white to black, bend the associated finger';
                task_label{ii} = 'fingermapping';
            elseif contains(task_label{ii}, 'gestures')        
                ieeg_json.TaskName = 'gestures';
                ieeg_json.TaskDescription = 'Make one of four hand gestures based on a learned visual cue';    
                ieeg_json.Instructions = 'When the visual cue appears, make the associated hand gesture';           
            end
            
        case 'tactile'
            % TACTILE TASKS
            if contains(task_label{ii}, 'temporalpattern') 
                ieeg_json.TaskName = 'temporalpattern';
                ieeg_json.TaskDescription = 'Tactile vibrations presented for variable durations and inter-stimulus intervals';
                ieeg_json.Instructions = 'Look at the fixation cross';
            end
            
        case 'none'
            % NO SPECIFIC TASKS - Skip
    end
    
    % Fix some problems in older stimulus files 
    if contains(task_label{ii}, 'hrf')   
        if length(find(stimData(ii).stimulus.trigSeq)) ~=32 && max(stimData(ii).stimulus.trigSeq) ~= 256
            % There are additional triggers that we should not use
            warning('[%s] Hrf stimulus file has incorrect number of triggers', mfilename)
            fprintf('[%s] Removing directly adjacent triggers from stimfile \n', mfilename);
            requestedTriggerInx = find(stimData(ii).stimulus.trigSeq);
            requestedTriggerInx_diff = [nan diff(requestedTriggerInx)];
            stimData(ii).stimulus.trigSeq(requestedTriggerInx(requestedTriggerInx_diff == 1)) = 0;
        end
    end
    if contains(task_label{ii}, 'prf')
        if length(find(stimData(ii).stimulus.trigSeq)) ~=224 && max(stimData(ii).stimulus.trigSeq) ~= 256
            % Old prf stim files had too many triggers + missing triggers
            warning('[%s] Prf stimulus file has incorrect number of triggers', mfilename)
            [newTrigSeq] = bidsconvert_fixprftrigseq(stimData(ii).stimulus);
            stimData(ii).stimulus.trigSeq = newTrigSeq;
            % Overwrite trial names to prevent concatenation problems in
            % preprocessing
            stimData(ii).stimulus.tsv.trial_name = [repmat('PRF', [height(stimData(ii).stimulus.tsv),1]) num2str(stimData(ii).stimulus.tsv.trial_name)];
        end
    end
    if contains(task_label{ii}, 'temporalpattern')
        if length(find(stimData(ii).stimulus.trigSeq)) ~= 144
            % There are additional triggers that we should not use
            warning('[%s] Temporal pattern stimulus file has incorrect number of triggers', mfilename)
            fprintf('[%s] Removing the triggers of block onsets/offsets from stimfile \n', mfilename);
            newTrigSeq = stimData(ii).stimulus.trigSeq;
            newTrigSeq(1) = 0;
            newTrigSeq(end) = 0;
            stimData(ii).stimulus.trigSeq = newTrigSeq;
        end
    end

    % Get the onsets in the data recordings for this run
    t = ((0:hdr.nSamples-1)/hdr.Fs); 
    num_triggers = length(find(stimData(ii).stimulus.trigSeq));
    num_triggers_total = num_triggers_total + num_triggers;
    onsets = trigger_onsets(stim0:num_triggers_total); 
    [~,onset_indices] = intersect(t, onsets);
   
    % Check if there are onset/offset triggers
    if (isequal(sensoryDomain, 'visual') && max(stimData(ii).stimulus.trigSeq) == 256) || isequal(sensoryDomain, 'motor') 
        % THESE DATA HAVE ONSET/OFFSET TRIGGERS
        hasOnOffTriggers = 1;
        % First trigger is task onset trigger:
        run_start_inx = onset_indices(1);
        % Last trigger is task offset trigger:
        run_stop_inx = onset_indices(end);
        % Stimulus onsets are those in between
        onset_indices = onset_indices(2:end-1);
    elseif isequal(sensoryDomain, 'tactile')
        % THESE DATA DO NOT HAVE ONSET/OFFSET TRIGGERS
        hasOnOffTriggers = 0;
        % Segment 3 seconds before first onset
        run_start_inx = max(onset_indices(1)-3*hdr.Fs,1);
        % Segment 3 seconds after last onset
        run_stop_inx = min(onset_indices(end)+3*hdr.Fs,size(data,2)); 
    elseif isequal(sensoryDomain, 'none')
        % THESE DATA DO NOT HAVE ONSET/OFFSET TRIGGERS
        hasOnOffTriggers = 0;
        % Segment from first onset
        run_start_inx = onset_indices(1);
        % Segment to end
        run_stop_inx = size(data,2)+1;
    else
        % THESE DATA DO NOT HAVE ONSET/OFFSET TRIGGERS
        hasOnOffTriggers = 0;
        % Segment 3 seconds before first onset
        run_start_inx = max(onset_indices(1)-3*hdr.Fs,1);
        % Segment 3 seconds after last onset
        run_stop_inx = min(onset_indices(end)+3*hdr.Fs,size(data,2));             
    end
    
    % Clip data from run using task onset and offset triggers
    data_thisrun = data(:,run_start_inx:run_stop_inx-1); % subtract one sample to make length of 'between run' data exactly 6 seconds
    hdr_thisrun = hdr;
    hdr_thisrun.nSamples = size(data_thisrun,2);
    
    if segmentOnFlips
        % Determine trial onset based on the flips:
        % Get the fliptimes for the requested triggers
        fprintf('Flip array length: %d | TrigSeq length: %d\n', ...
            numel(stimData(ii).response.flip), numel(stimData(ii).stimulus.trigSeq));

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
    
    % Collect info for events.tsv file
    events_table = stimData(ii).stimulus.tsv;
    
    % Overwrite onset with onsets of triggers
    events_table.event_sample = event_sample;
    events_table.onset        = strtrim(cellstr(num2str(events_table.event_sample/hdr.Fs,['%.' num2str(nDecimals) 'f'])));
    
    % Update a number of other fields in events table:
    
    % Some formatting to facilitate concatenation of events files in analysis:
    if ~isfield(summary(events_table), 'ISI'); events_table.ISI = zeros(height(events_table),1); end
    events_table.duration     = strtrim(cellstr(num2str(events_table.duration,['%.' num2str(nDecimals) 'f']))); 
    events_table.ISI          = strtrim(cellstr(num2str(events_table.ISI,['%.' num2str(nDecimals) 'f'])));    
    % Add a task column to the events_table 
    events_table.task_name   = repmat(task_label{ii}, height(events_table), 1);
    
    % Collect info for json_ieeg file
    ieeg_json.SamplingFrequency = hdr_thisrun.Fs;
    ieeg_json.RecordingDuration = round(hdr_thisrun.nSamples/hdr_thisrun.Fs,nDecimals);
  
    %%  Write data, channel, json and events file:    
    
    % Generate a filename
    fname = sprintf('sub-%s_ses-%s_task-%s_acq-%s_run-%s', ...
            sub_label, ses_label, task_label{ii}, acq_label, run_label{ii});
        
    fprintf('[%s] Writing bids files for %s \n', mfilename, fname);
    
    % Write data files:    
    data_fname = fullfile(dataWriteDir, sprintf('%s_ieeg', fname));
    ft_write_data(data_fname, data_thisrun, 'header', hdr_thisrun, 'dataformat', 'brainvision_eeg');
 
    % Write json_ieeg file:
    jsonfile_fname = fullfile(dataWriteDir, sprintf('%s_ieeg.json', fname));    
    jsonwrite(jsonfile_fname,ieeg_json,json_options)
    
    % Write channels.tsv file:
    channels_fname = fullfile(dataWriteDir, sprintf('%s_channels.tsv', fname));    
    writetable(channel_table,channels_fname,'FileType','text','Delimiter','\t');
    
	% Write stimulus file:
    stimfile_thisrun = stimData(ii);
    stimfile_fname = fullfile(stimWriteDir, sprintf('%s.mat', fname));
    save(stimfile_fname,'-struct', 'stimfile_thisrun', '-v7.3')
    
    % Update events table to point to the newly saved stimulus file:
    events_table.stim_file = repmat([fname '.mat'], height(events_table), 1);
    
    % Write events.tsv file: 
    events_fname = fullfile(dataWriteDir, sprintf('%s_events.tsv', fname));
    writetable(events_table, events_fname, 'FileType','text', 'Delimiter', '\t')
    
    % Collect filenames to be used for scants.tsv in
    % bidsconvert_writesessionfiles.m
    fname = fullfile('ieeg',sprintf('%s_ieeg.eeg', fname));
	dataFileNames{ii} = fname;

end

% CHECK: Do number of triggers derived from EDF data file match the number
% of trials from the stimulus files?
assert(isequal(length(trigger_onsets), num_triggers_total))
fprintf('[%s] BIDS conversion done\n', mfilename)

end