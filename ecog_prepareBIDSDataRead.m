function [data, events, chans] = ecog_readBIDSData(dataDir, sub_label, ses_label)

% Specify fieldtrip configuraton for reading in the files
cfg            = [];
cfg.continuous = 'yes';
cfg.channel    = 'all';

cfg.dataset = [];
secondsToAdd = 0;

% Identify the number of runs 
dataFiles   = dir(fullfile(dataDir,sprintf('sub-%s_ses-%s_task-*.eeg',sub_label, ses_label)));
eventFiles  = dir(fullfile(dataDir,sprintf('sub-%s_ses-%s_task-*events.tsv',sub_label, ses_label)));
chanFiles   = dir(fullfile(dataDir,sprintf('sub-%s_ses-%s_task-*channels.tsv',sub_label, ses_label)));

for iRun = 1:length(dataFiles)
    
    % Generate a list of datafiles to read in with ft_preprocessing below
    cfg.dataset = [cfg.dataset {fullfile(dataDir,dataFiles(iRun).name)}];
      
    % Read in events for this run
    eventsTable = readtable(fullfile(dataDir,eventFiles(iRun).name), 'FileType', 'text');
    
    if iRun == 1 
        events = eventsTable;
    else
        eventsTable.onset = eventsTable.onset+secondsToAdd; % add length of PREVIOUS run in seconds
        events = [events; eventsTable];
    end
    
    % Read in hdr for this run; use to update event onsets for next run
    hdr = ft_read_header(fullfile(dataDir,dataFiles(iRun).name));
    runLengthInSeconds = round(hdr.nSamples/hdr.Fs, nDecimals);
    secondsToAdd = secondsToAdd + runLengthInSeconds;
    fprintf('Length of run %d is %s seconds, cumulative length is %s seconds \n', iRun, num2str(runLengthInSeconds), num2str(secondsToAdd));
end

% Read in one of the channel files
chans   = readtable(fullfile(dataDir,chanFiles(1).name), 'FileType', 'text');


