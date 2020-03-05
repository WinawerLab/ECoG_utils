function [data, channels, events] = bidsEcogGetPreprocData(dataPath, subject, sessions, tasks, runnums, description)
% Reads in timeseries data, channels and events from a BIDS directory with
% ECoG data. 
%   
% [data, channels, events, srate] = bidsEcogGetPreProcData(dataPath, subject, [sessions], [tasks], [runnums])
%
% Inputs
%   dataPath:   path to folder containing preprocessed data
%   dataStr:    text string to specify filename for data
%   tasks:      BIDS tasks, in cell array
%   runnums:    cell array of runnumbers, equal in length to tasks
%
% Output
%   data:       the time-series data for each run with dimensions
%                X x Y x Z x time
%   info:       nifti header for each run
%   fullFile:   The (nifti) time-series data including header information
%
% Example:

[sessions] = bidsSpecifySessions(dataPath, subject, sessions);
if ~iscell(sessions), sessions = {sessions}; end

% Initialize
allData = []; 
allEvents = [];
samplesToAdd = 0;
secondsToAdd = 0;

if isempty(runnums), getRunsPerSession = 1; else, getRunsPerSession = 0; end
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

                % Run various checks on events table:

                % Check if there are trial_names, if not, add them
                if ~isfield(summary(events), 'trial_name')
                    fprintf('[%s] Events lack trial names - adding them now .\n',mfilename);
                    events = bair_addTrialNamesToEventsTable(events);
                end
                % Check if there are trial_names, if not, add them
                if ~isfield(summary(events), 'stim_file')
                    fprintf('[%s] Events lack stim files - adding them now .\n',mfilename);
                    events.stim_file = repmat({'n/a'}, [height(events) 1]);
                end

                % Check if there are sample indices and ISIs:
                if ~isfield(summary(events),'event_sample'); events.event_sample = round(events.onset*hdr.Fs);end
                if ~isfield(summary(events), 'ISI');events.ISI = zeros(height(events),1);end

                % Add task, session and run indices to the events file
                if ~isfield(summary(events), 'task_name'), events.task_name = repmat({task}, [height(events),1]); end
                events.session_name = repmat({session}, [height(events),1]);
                events.run_name = repmat({runnum}, [height(events),1]);

                % Concatenate data and events; update onsets 
                if runCount == 1 
                    allEvents = events; 
                    allData = data; 
                end
                if runCount > 1
                    events.event_sample = events.event_sample + samplesToAdd;
                    events.onset = events.onset + secondsToAdd;
                    allEvents = [allEvents; events];                    
                    allData = cat(2,allData,data);
                end

                % Use runlengths to update next run
                runLengthInSamples = hdr.nSamples;
                runLengthInSeconds = hdr.nSamples/hdr.Fs;
                samplesToAdd = samplesToAdd + runLengthInSamples;
                secondsToAdd = secondsToAdd + runLengthInSeconds;
                fprintf('[%s] Length of run %s is %s seconds, cumulative length of data is %s seconds \n', mfilename, runnum, num2str(runLengthInSeconds), num2str(samplesToAdd/hdr.Fs));          
            end
        end
        % Check channel statuses
        if ~exist('previousChannels', 'var') 
            previousChannels = channels; 
        else
            if isfield(summary(channels), 'status')
                if ~isequal(channels.status, previousChannels.status)  
                     % This means one or more electrodes are bad in one session
                     % but not another. For now, label the electrode as good
                     % across all sessions:
                     channels.status(~strcmp(channels.status, previousChannels.status)) = {'good'};
                     % and assume bad channels/epochs will removed at later
                     % stages.
                     % TO DO: find some smarter way to deal with this problem
                end
            end
        end
    end
end

% Remove non-data channels.
chan_inx = contains(lower(channels.type), {'ecog', 'seeg'}); 
channels = channels(chan_inx,:);
data = allData(chan_inx,:);
events = allEvents;

if exist('channels', 'var') 
    channels.subject_name = repmat({subject},[height(channels) 1]);
else
    channels = [];
end

end