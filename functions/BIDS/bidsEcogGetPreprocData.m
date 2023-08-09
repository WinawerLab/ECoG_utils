function [data, channels, events, srate, nsamples] = bidsEcogGetPreprocData(dataPath, subject, sessions, tasks, runnums, description, srate)
% Reads in timeseries data, channels and events from a BIDS directory with
% preprocessed ECoG data. Data from specified tasks and runs will be
% concatenated, and the event onsets will be updated accordingly.
%   
% [data, channels, events, srate, nsamples] = bidsEcogGetPreProcData(dataPath, ...
%           subject, sessions, tasks, runnums, [description], [srate])
%
% Inputs
%     dataPath:         path where the BIDS project lies (string) 
%     subject:          BIDS subject name (string, all lower case) 
%     sessions:         [optional] BIDS session name (cell of strings, all lower
%                       case). default: all sessions for this subject.
%     tasks:            [optional] BIDS task name (cell of strings, all
%                       lower case). default: all tasks for this subject.
%     runnums:          [optional] BIDS run number (cell of strings, all
%                       lower case). default: all runs for this subject.
%     description:      [optional] String for the BIDS 'desc-' label in the
%                       bidsname. default: all desc- labels in dataPth.
%     srate:            [optional] Integer indicating a desired sampling
%                       rate in Hz; data that does not match sampling rate
%                       will be resampled to this rate. Useful for when
%                       different runs are recorded with different srates.
%
% Output
%     data:             A 2D matrix with ECoG data (channels x time),
%                       concatenated across tasks and runs
%     channels:         A table with channel information (concatenated
%                       across tasks and runs)
%     events:           A table with events information (concatenated
%                       across tasks and runs, with updated onsets)
%     srate:        	Sampling rate (native when the no srate input is
%                       provided, should match input srate if provided).
%     nsamples:         Number of samples for loaded files
%
% Example:
%  dataPath = fullfile(bidsRootPath, 'derivatives', 'ECoGCAR');
%  subject = 'p08'; 
%  tasks = {'spatialpattern', 'temporalpattern'};
%  [data, channels, events, srate] = bidsEcogGetPreprocData(dataPath, subject, [], tasks);
      
% <dataPath>
if ~exist('dataPath', 'var') || isempty(dataPath)
    error('dataPath not defined');
end  

% <subject>
if ~exist('subject', 'var') || isempty(subject)
    error('subject not defined');
end

% <session>
if ~exist('sessions', 'var') || isempty(sessions)
	[sessions] = bidsSpecifySessions(dataPath, subject);
end

% <task>
if ~exist('tasks', 'var'), tasks = []; end

% <runnum>
if ~exist('runnums', 'var') || isempty(runnums)
	getRunsPerSession = 1; 
else
    getRunsPerSession = 0;
end

% <description>
if ~exist('description', 'var'), description = []; end

% <srate>
if ~exist('srate', 'var'), srate = []; end

if ~iscell(sessions), sessions = {sessions}; end

% Initialize
allData = []; 
allEvents = [];
samplesToAdd = 0;
secondsToAdd = 0;
nsamples = [];
runCount = 0;

for ii = 1:length(sessions)

    % check how many runs there are for each task
    if getRunsPerSession, runnums = []; end
    [session, tasks, runnums] = bidsSpecifyData(dataPath, subject, sessions{ii}, tasks, runnums);
     
    for jj = 1:length(tasks)

        task = tasks{jj};
        
        % check if we have any runs this task in this session; 
        % if not, continue to next session
    
        if isempty(runnums{jj})
            
            fprintf('[%s] No runs for task %s in session %s \n',mfilename, task, session);
            continue
            
        else
            fprintf('[%s] Concatenating runs for session %s, task %s...\n',mfilename, session, task);

            for kk = 1:length(runnums{jj})

                runnum = runnums{jj}{kk};

                [data, channels, events, ~, hdr] = bidsEcogReadFiles(dataPath, subject, session, task, runnum, description);

                runCount = runCount + 1;
            
                % Check sampling frequency
                if ~isempty(srate)
                    Fs = hdr.Fs;
                    if ~isequal(Fs,srate)
                        fprintf('[%s] Sample rate of %d Hz does not match requested sample rate of %d Hz. Resampling... \n',mfilename, hdr.Fs, srate);
                        data = resample(data', srate, hdr.Fs)';
                        channels.sampling_frequency(:) = srate;
                    end
                else
                    if exist('previousFs', 'var') || ~isempty(previousFs)
                        assert(hdr.Fs == previousFs,'Failed to concatenate sessions because the sampling frequencies were different. Please specify a desired sampling rate using the input argument srate.');
                    end
                end
                previousFs = hdr.Fs;

                % Check channel statuses
                if exist('previousChannels', 'var')
                    if isfield(summary(channels), 'status')
                        if ~isequal(channels.status, previousChannels.status)
                            % This means one or more electrodes are bad in
                            % one session but not another. For now, label
                            % the electrode as good across all sessions,
                            % and assume bad epochs will removed later.
                            channels.status(~strcmp(channels.status, previousChannels.status)) = {'good'};        
                        end
                    end
                end
                previousChannels = channels;

                % Run various checks on events table:

                % Check if there are trial_names, if not, add them
                if ~isfield(summary(events), 'trial_name')
                    if isnumeric(events.trial_type(1))
                        fprintf('[%s] Events lack trial names - adding them now .\n',mfilename);
                        events = bair_addTrialNamesToEventsTable(events);
                    else
                        events.trial_name = events.trial_type;
                    end
                end
                
                % Check for value field
                if isfield(summary(events), 'value')
                    if ~iscell(events.value)
                        events.value = num2cell(events.value);
                    end
                end
                
                % Check if there are stim_file_names, if not, add them
                if ~isfield(summary(events), 'stim_file')
                    fprintf('[%s] Events lack stim files - adding them now .\n',mfilename);
                    events.stim_file = repmat({'n/a'}, [height(events) 1]);
                end

                % Check if there are ISIs:
                if ~isfield(summary(events), 'ISI');events.ISI = zeros(height(events),1);end

                % Add task, session and run indices to the events file
                if ~isfield(summary(events), 'task_name'), events.task_name = repmat({task}, [height(events),1]); end
                events.session_name = repmat({session}, [height(events),1]);
                events.run_name = repmat({runnum}, [height(events),1]);

                % Concatenate data and events; update onsets 
                if runCount > 1
                    events.onset = events.onset + secondsToAdd;
                    allEvents = [allEvents; events];                    
                    allData = cat(2,allData,data);
                else
                    allEvents = events; 
                    allData = data; 
                end
                nsamples = cat(1,nsamples,length(data));

                % Use runlengths to update next run
                runLengthInSamples = hdr.nSamples;
                runLengthInSeconds = hdr.nSamples/hdr.Fs;
                if isempty(runLengthInSamples), runLengthInSamples = 0; end
                if isempty(runLengthInSeconds), runLengthInSeconds = 0; end
                samplesToAdd = samplesToAdd + runLengthInSamples;
                secondsToAdd = secondsToAdd + runLengthInSeconds;
                fprintf('[%s] Length of run %s is %0.1f seconds, cumulative length of data is %0.1f seconds \n', mfilename, runnum, runLengthInSeconds, samplesToAdd./hdr.Fs);          
            end
        end        
    end
end

% Remove non-data channels.
if exist('channels', 'var') 
    chan_inx = contains(lower(channels.type), {'ecog', 'seeg'}); 
    channels = channels(chan_inx,:);
    channels.subject_name = repmat({subject},[height(channels) 1]);
    data = allData(chan_inx,:);
    events = allEvents;
    if isempty(srate), srate = hdr.Fs; end
else
    channels = [];
    data     = [];
    events   = []; 
    srate    = [];
    nsamples = [];
end

end