
function [data, channels, events, json] = bidsEcogGetPreprocData(dataPath, subject, sessions, tasks, runnums, description)

% get timeseries, channels, events, concatenated for subset of tasks
% should get timeseries, events and channel info for a subset of
% tasks and runnums(optional) for a given data folder.
    
% [data, channels, events] = bidsEcogGetPreProcData(dataPath, subject, [sessions], [tasks], [runnums])

[sessions] = bidsSpecifySessions(dataPath, subject, sessions);

% Initialize
allData = []; 
allEvents = [];
samplesToAdd = 0;
secondsToAdd = 0;
    
for ii = 1:length(sessions)

    % bidsSpecifyData
    [session, tasks, runnums] = bidsSpecifyData(dataPath, subject, sessions{ii}, tasks, runnums);
     
    for jj = 1:length(tasks)
        
        task = tasks{jj};
        
        fprintf('[%s] Concatenating runs for session %s, task %s...\n',mfilename, session, task);
        
        for kk = 1:length(runnums{jj})
            
            runnum = runnums{jj}{kk};

            [data, channels, events, json, hdr] = bidsEcogReadFiles(dataPath, subject, session, task, runnum, description);
            
            % Run various checks:
            
            % Check if there are trial_names, if not, add them
            if ~isfield(summary(events), 'trial_name')
                fprintf('[%s] Event lack trial names - adding them now .\n',mfilename);
                events = bair_addTrialNamesToEventsTable(events);
            end
            
            if ~iscell(events.stim_file_index)
                events.stim_file_index = num2cell(events.stim_file_index);
            end
 
            if ~isfield(events, 'ISI')
                events.ISI = zeros(height(events),1);
            end
            
            % Add task and run indices to the events file
            events.run_name = repmat({runnum}, [height(events),1]);
            events.session_name = repmat({session}, [height(events),1]);
            
            % Concatenate tasks and runs - add task and run indices to events 

            % Concatenate data and events; update onsets 
            if jj == 1, allEvents = events; end
            if jj > 1
                if isfield(summary(events),'event_sample')
                    events.event_sample = events.event_sample + samplesToAdd;
                    events.onset = events.event_sample/hdr.Fs;
                else
                    events.onset = events.onset + secondsToAdd;
                end
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
    
    % When concatenating across sessions, check channel statuses
    if ii == 1, previousChannels = channels; end
    if ii > 1
        if isfield(summary(channels), 'status')
            assert(isequal(channels.status, previousChannels.status))
%                 % todo: if this assertion fails, meaning that a visual
%                 % electrode is bad in one session but not another, include
%                 % only the good session, e.g. put one set of trials to nans
        end
    end
end
data = allData;
events = allEvents;