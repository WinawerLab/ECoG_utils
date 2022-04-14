function [data, channels, events, ieeg_json, hdr] = bidsEcogReadFiles(dataPath, subject, session, task, runnum, description)
% Reads in ECoG data files from dataPath with BIDS formatted data. 
%
% [data, channels, events, ieeg_json, hdr] = bidsEcogReadFiles(dataPath, ...
%     subject, session, task, runnum, [description])
%
% Input
%     dataPath:         path where the BIDS project lies (string)
%     subject:          BIDS subject name (string, all lower case)
%     session:          BIDS session name (string, all lower case)
%     task:             BIDS task name (string, all lower case)                           
%     runnum:           BIDS run number (string, all lower case)
%     description:      String for the BIDS 'desc-' label in the bidsname
%                       [optional], if not provided or empty, 'desc-' is 
%                       not added to bidsname.
%
% Output
%     data:             A 2D matrix with ECoG data (channels x time)
%     channels:         A table with channel information
%     events:           A table with events information
%     ieeg_json:        A struct with the ECoG meta data
%     hdr:              The data header as provided by Fieldtrip
%
% Note
% This function uses FieldTrips ft_read_data.m and ft_read_hdr.m to read in
% the ECOG data. Additional expected meta-data files are a channels file
% (<bidsfilename>_channels.tsv), an events file (<bidsfilename>_events.tsv)
% and a json sidecar (<bidsfilename>_ieeg.json). 
%
% Example
%   projectDir = '/Volumes/server/Projects/BAIR/Data/BIDS/visual';
%   dataPath = fullfile(projectDir, 'derivatives', 'ECoGBroadband');
%   subject = 'p10';
%   session = 'nyuecog01'
%   task = 'prf';
%   runnum = '01';
%   description = 'broadband';
%  [data, channels, events, ieeg_json, hdr] = bidsEcogReadFiles(dataPath, ...
%       subject, session, task, runnum, description);
%
% See also bidsEcogWriteFiles.m bidsSpecifyData.m bidsSpecifySessions.m
% See also ft_read_data.m ft_read_hdr.m
%
% IG, BAIR 2019

% <dataPath>
if ~exist('dataPath', 'var') || isempty(dataPath)
    error('dataPath not defined');
end  

% <subject>
if ~exist('subject', 'var') || isempty(subject)
    error('subject not defined');
end

% <session>
if ~exist('session', 'var') || isempty(session)
	error('session not defined');
end

% <task>
if ~exist('task', 'var') || isempty(task)
	error('task not defined');
end

% <runnum>
if ~exist('runnum', 'var') || isempty(runnum)
	error('runnum not defined');
end

% <description>
if ~exist('description', 'var'), description = []; end

if ~isempty(description)
    fname_in = sprintf('sub-%s_ses-%s_task-%s*_run-%s_desc-%s*', subject, session, task, runnum, description);
else
    fname_in = sprintf('sub-%s_ses-%s_task-%s*_run-%s*', subject, session, task, runnum);
end

sessionDir = fullfile(dataPath, sprintf('sub-%s', subject), sprintf('ses-%s', session), 'ieeg');
fprintf('[%s] Reading from %s\n', mfilename, sessionDir);
fprintf('[%s] Reading %s\n', mfilename, fname_in);

% Read in the channels file
chanFile = dir(fullfile(sessionDir, sprintf('%s_channels.tsv', fname_in)));
if length(chanFile) < 1, error('did not find channels file %s in %s', fname_in, sessionDir); end
if length(chanFile) > 1, error('found multiple channels files %s in %s', fname_in, sessionDir); end
%fprintf('[%s] Reading in channels file: %s\n', mfilename, chanFile); 
fprintf('.'); 
channels = readtable(fullfile(chanFile.folder, chanFile.name), 'FileType', 'text');
%channels = tdfread(chanFile);
%channels = struct2table(channels);

% Read in the events file
eventsFile = dir(fullfile(sessionDir, sprintf('%s_events.tsv', fname_in)));
if length(eventsFile) < 1, error('did not find event file %s in %s', fname_in, sessionDir); end
if length(eventsFile) > 1, error('found multiple event files %s in %s', fname_in, sessionDir); end
%fprintf('[%s] Reading in events file: %s\n', mfilename, eventsFile); 
fprintf('.'); 
events = readtable(fullfile(eventsFile.folder, eventsFile.name), 'FileType', 'text');
%events = tdfread(eventsFile);
%events = struct2table(events);
 
% Read in the data file
dataFile = dir(fullfile(sessionDir, sprintf('%s_ieeg.eeg', fname_in)));
if length(dataFile) < 1, error('did not find data file %s in %s', fname_in, sessionDir); end
if length(dataFile) > 1, error('found multiple data files %s in %s', fname_in, sessionDir); end
%fprintf('[%s] Reading in data file: %s\n', mfilename, dataFile); 
fprintf('.'); 
hdr = ft_read_header(fullfile(dataFile.folder, dataFile.name));
data = ft_read_data(fullfile(dataFile.folder, dataFile.name));

% Read in the json file
jsonFile = dir(fullfile(sessionDir, sprintf('%s_ieeg.json', fname_in)));
if length(jsonFile) < 1, error('did not find json file %s in %s', fname_in, sessionDir); end
if length(jsonFile) > 1, error('found multiple json files %s in %s',fname_in, sessionDir); end
%fprintf('[%s] Reading in ieeg json file: %s\n', mfilename, jsonReadFile); 
fprintf('.\n'); 
ieeg_json = jsonread(fullfile(jsonFile.folder, jsonFile.name));
                                   
end

