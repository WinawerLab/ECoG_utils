
function [ts, channels, events] = bidsEcogGetPreProcData(dataPath, subject, sessions, tasks, runnums)

% get timeseries, channels, events, concatenated for subset of tasks
% should get timeseries, events and channel info for a subset of
% tasks and runnums(optional) for a given data folder.
    
% [data, channels, events] = bidsEcogGetPreProcData(dataPath, subject, [sessions], [tasks], [runnums])

[sessions] = bidsSpecifySessions(dataPath, subject, sessions);
        
for ii = 1:length(sessions)

    % bidsSpecifyData
    [session, tasks, runnum] = bidsSpecifyData(projectDir, subject, sessions{ii}, tasks);
    
    % get time series
    [ts, info] = bidsGetPreprocData(dataPath, dataStr, tasks, runnum);
    
    % get channels and events
    
    % do some checks, e.g. does channel table match the size of the ts
    
    % concatenate tasks and runs - add task and run indices to events 
    
end