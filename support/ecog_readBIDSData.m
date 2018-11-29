function [data, events, chans] = ecog_readBIDSData(dataDir, sub_label, ses_label)

% Specify fieldtrip configuraton for reading in the files
cfg            = [];
cfg.continuous = 'yes';
cfg.channel    = 'all';

% Identify the number of runs 
dataFiles   = dir(fullfile(dataDir,sprintf('sub-%s_ses-%s_task-*.eeg',sub_label, ses_label)));
eventFiles  = dir(fullfile(dataDir,sprintf('sub-%s_ses-%s_task-*events.tsv',sub_label, ses_label)));
chanFiles   = dir(fullfile(dataDir,sprintf('sub-%s_ses-%s_task-*channels.tsv',sub_label, ses_label)));

% Pre-allocate list of runs to be concatenated 
cfg.dataset = [];
samplesToAdd = 0;

fprintf('[%s] Concatenating runs...\n',mfilename);

% Loop across all the runs
for iRun = 1:length(dataFiles)
    
    % Generate a list of datafiles to read in with ft_preprocessing below
    cfg.dataset = [cfg.dataset {fullfile(dataDir,dataFiles(iRun).name)}];
      
    % Read in events for this run
    eventsTable = readtable(fullfile(dataDir,eventFiles(iRun).name), 'FileType', 'text');
    
    % Concatenate events across runs
    if iRun == 1 
        events = eventsTable;
    else
        % Add length of PREVIOUS run in seconds to onsets
        eventsTable.event_sample = eventsTable.event_sample + samplesToAdd;
        eventsTable.onset = eventsTable.event_sample/hdr.Fs;
        events = [events; eventsTable];
    end
    
    % Read in hdr for this run; use to update event onsets for next run
    hdr = ft_read_header(fullfile(dataDir,dataFiles(iRun).name));
    runLengthInSamples = hdr.nSamples;
    runLengthInSeconds = hdr.nSamples/hdr.Fs;
    samplesToAdd = samplesToAdd + runLengthInSamples;
    fprintf('[%s] Length of run %d is %s seconds, cumulative length of data is %s seconds \n', mfilename,iRun, num2str(runLengthInSeconds), num2str(samplesToAdd/hdr.Fs));
end

% Read in the data files (all runs)
data    = ft_preprocessing(cfg);

% Read in one of the channel files
chans   = readtable(fullfile(dataDir,chanFiles(1).name), 'FileType', 'text');

% Check whether the channel info in channels.tsv matches the data
assert(size(chans,1) == length(data.label));

% Replace data.label with names from channels.tsv for better readability
data.label = chans.name;
data.hdr.label = chans.name;
%fprintf('done \n');

end
