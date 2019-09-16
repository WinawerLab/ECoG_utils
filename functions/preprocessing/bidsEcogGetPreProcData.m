
function [ts, channels, events] = bidsEcogGetPreprocData(dataPath, subject, sessions, tasks, runnums, description)

% get timeseries, channels, events, concatenated for subset of tasks
% should get timeseries, events and channel info for a subset of
% tasks and runnums(optional) for a given data folder.
    
% [data, channels, events] = bidsEcogGetPreProcData(dataPath, subject, [sessions], [tasks], [runnums])

[sessions] = bidsSpecifySessions(dataPath, subject, sessions);
        
for ii = 1:length(sessions)

    % bidsSpecifyData
    [session, tasks, runnums] = bidsSpecifyData(dataPath, subject, sessions{ii}, tasks, runnums);
    
    for ii = 1:length(tasks)
        for jj = 1:length(runnums{ii})
            
            
            [data, channels, events] = bidsEcogReadFiles(dataPath, subject, session, tasks{ii}, runnums{ii}{jj}, description);
            

        end
    end

    
     
%     % Pre-allocate list of runs to be concatenated 
%     cfg.dataset = [];
%     samplesToAdd = 0;
%     secondsToAdd = 0;
% 
%     fprintf('[%s] Concatenating runs...\n',mfilename);
% 
%     % Loop across all the runs
%     for iRun = 1:length(dataFiles)
% 
%         % Generate a list of datafiles to read in with ft_preprocessing below
%         cfg.dataset = [cfg.dataset {fullfile(dataDir,dataFiles(iRun).name)}];
% 
%         % Read in events for this run
%         eventsTable = readtable(fullfile(dataDir,eventFiles(iRun).name), 'FileType', 'text');
%         if isfield(summary(eventsTable), 'trial_name')
%             if ~iscell(eventsTable.trial_name)
%                 if contains(eventFiles(iRun).name, 'prf')
%                     eventsTable.trial_name = repmat({'PRF'}, [height(eventsTable) 1]);
%                 elseif contains(eventFiles(iRun).name, 'hrf')
%                     eventsTable.trial_name = repmat({'HRF'}, [height(eventsTable) 1]);
%                 end
%             end
%         end
% 
%         if ~iscell(eventsTable.stim_file_index)
%              eventsTable.stim_file_index = num2cell(eventsTable.stim_file_index);
%         end
% 
%         if ~isfield(eventsTable, 'ISI')
%             eventsTable.ISI = zeros(height(eventsTable),1);
%         end
% 
%         % Concatenate events across runs
%         if iRun == 1 
%             events = eventsTable;
%         else
%             % Add length of PREVIOUS run in seconds to onsets
%             %if isfield(summary(eventsTable),'event_sample')
%             if max(contains(eventsTable.Properties.VariableNames, 'event_sample'))>0
%                 eventsTable.event_sample = eventsTable.event_sample + samplesToAdd;
%                 eventsTable.onset = eventsTable.event_sample/hdr.Fs;
%             else
%                 eventsTable.onset = eventsTable.onset + secondsToAdd;
%             end
%             events = [events; eventsTable];
%         end
% 
%         % Read in hdr for this run; use to update event onsets for next run
%         hdr = ft_read_header(fullfile(dataDir,dataFiles(iRun).name));
%         runLengthInSamples = hdr.nSamples;
%         runLengthInSeconds = hdr.nSamples/hdr.Fs;
%         samplesToAdd = samplesToAdd + runLengthInSamples;
%         secondsToAdd = secondsToAdd + runLengthInSeconds;
%         fprintf('[%s] Length of run %d is %s seconds, cumulative length of data is %s seconds \n', mfilename,iRun, num2str(runLengthInSeconds), num2str(samplesToAdd/hdr.Fs));
%     end
% 
%     fprintf('[%s] Reading in the data using FieldTrip:\n',mfilename);
%     % Read in the data files (all runs)
%     data    = ft_preprocessing(cfg);
% 
%     % Read in one of the channel files
%     chans   = readtable(fullfile(dataDir,chanFiles(1).name), 'FileType', 'text');
% 
%     % Check whether the channel info in channels.tsv matches the data
%     assert(size(chans,1) == length(data.label));
% 
%     % Replace data.label with names from channels.tsv for better readability
%     data.label = chans.name;
%     data.hdr.label = chans.name;
%     %fprintf('done \n');
    % do some checks, e.g. does channel table match the size of the ts
    
    % concatenate tasks and runs - add task and run indices to events 
    
end