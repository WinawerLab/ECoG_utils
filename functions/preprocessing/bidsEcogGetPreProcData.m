
function [ts, channels, events] = bidsEcogGetPreProcData(dataPath, subject, sessions, tasks, runnums)

% get timeseries, channels, events, concatenated for subset of tasks

% [data, channels, events] = bidsEcogGetPreProcData(dataPath, subject, [sessions], [tasks], [runnums])


for j = 1 : length(ses_label{k})




    % bidsSpecifyData

    % should get timeseries, events and channel info for a subset of
    % tasks and runnums(optional) for a given preprocessed folder.
    % runs and tasks(?) should be concatenated?

    % use? [data, info] = bidsGetPreprocData(dataPath, dataStr, tasks, runnums);



end