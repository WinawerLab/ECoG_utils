function [data, channels, events, ieeg_json, hdr] = bidsEcogReadFiles(dataPath, subject, session, task, runnum, description)

% TO DO add description, what to with less inputs

% <outputFolder>
if ~exist('description', 'var')
    description = [];
end

fprintf('[%s] Task = %s, Run = %s \n', mfilename, task, runnum);

if ~isempty(description)
    fname_in = sprintf('sub-%s_ses-%s_task-%s_run-%s_desc-%s', subject, session, task, runnum, description);
else
    fname_in = sprintf('sub-%s_ses-%s_task-%s_run-%s', subject, session, task, runnum);
end

sessionDir = fullfile(dataPath, sprintf('sub-%s', subject), sprintf('ses-%s', session));

% Read in the channels file
chanFile = fullfile(sessionDir, 'ieeg', sprintf('%s_channels.tsv', fname_in));
if ~exist(chanFile, 'file'), error('channels file not found: %s', chanFile); end
fprintf('[%s] Reading in channels file: %s\n', mfilename, chanFile); 
channels = readtable(chanFile, 'FileType', 'text');
%channels = tdfread(chanFile);
%channels = struct2table(channels);

% Read in the events file
eventsFile = fullfile(sessionDir, 'ieeg', sprintf('%s_events.tsv', fname_in));
if ~exist(eventsFile, 'file'), error('events file not found: %s', eventsFile); end
fprintf('[%s] Reading in events file: %s\n', mfilename, eventsFile); 
events = readtable(eventsFile, 'FileType', 'text');
%events = tdfread(eventsFile);
%events = struct2table(events);
 
% Read in the data file
dataFile = fullfile(sessionDir, 'ieeg', sprintf('%s_ieeg.eeg', fname_in));
if ~exist(dataFile, 'file'), error('data file not found: %s', dataFile); end
fprintf('[%s] Reading in data file: %s\n', mfilename, dataFile); 
hdr = ft_read_header(dataFile);
data = ft_read_data(dataFile);

% Read in the json file
jsonReadFile = fullfile(sessionDir, 'ieeg', sprintf('%s_ieeg.json', fname_in));
if ~exist(jsonReadFile, 'file'), error('data file not found: %s', dataReadFile); end
fprintf('[%s] Reading in ieeg json file: %s\n', mfilename, jsonReadFile); 
ieeg_json = jsonread(jsonReadFile);
           
                           
end

