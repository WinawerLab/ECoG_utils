function [data, channels, events, ieeg_json, hdr] = bidsEcogReadFiles(dataPath, subject, session, task, runnum, description)
% Reads in ECoG data files from dataPath with bidsfilename formatted as:
% sub-<subject>_ses-<session>_task-<task>_run-<runnum>_desc-<description>.
%
% [data, channels, events, ieeg_json, hdr] = bidsEcogReadFiles(dataPath, ...
% subject, session, task, runnum, [description])
%
%
% Input
%     dataPath:         path where the BIDS project lies (string)
%     subject:          BIDS subject name (string, all lower case)
%     session:          BIDS session name (string, all lower case)
%     task:             BIDS task name (string, all lower case)                           
%     runnum:           BIDS run number (string, all lower case)
%     description:      String for the BIDS 'desc-' label in the bidsname
%                       [optional], if not providedor empty, 'desc-' is not
%                       added to bidsname.
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
% the ECOG data. In addition, Expected meta-data files are a channels file
% (<bidsfilename>_channels.tsv), an events file (<bidsfilename>_events.tsv)
% and a json file (<bidsfilename>_ieeg.json). Missing files result in
% errors.
%
% Example
% projectDir = '/Volumes/server/Projects/BAIR/Data/BIDS/visual';
% dataPath = fullfile(projectDir, 'derivatives', 'ECoGBroadband');
% subject = 'som726';
% session = 'nyuecog01'
% task = 'prf';
% runnum = '01';
% description = 'broadband';
% [data, channels, events, ieeg_json, hdr] = bidsEcogReadFiles(dataPath, ...
% subject, session, task, runnum, description);
%
% See also bidsEcogWriteFiles.m bidsSpecifyData.m bidsSpecifySessions.m
% See also ft_read_data.m ft_read_hdr.m


if ~exist('description', 'var')
    description = [];
end

if ~isempty(description)
    fname_in = sprintf('sub-%s_ses-%s_task-%s_run-%s_desc-%s', subject, session, task, runnum, description);
else
    fname_in = sprintf('sub-%s_ses-%s_task-%s_run-%s', subject, session, task, runnum);
end

sessionDir = fullfile(dataPath, sprintf('sub-%s', subject), sprintf('ses-%s', session));
fprintf('[%s] Reading from %s\n', mfilename, sessionDir);
fprintf('[%s] Reading %s', mfilename, fname_in);

% Read in the channels file
chanFile = fullfile(sessionDir, 'ieeg', sprintf('%s_channels.tsv', fname_in));
if ~exist(chanFile, 'file'), error('channels file not found: %s', chanFile); end
%fprintf('[%s] Reading in channels file: %s\n', mfilename, chanFile); 
fprintf('.'); 
channels = readtable(chanFile, 'FileType', 'text');
%channels = tdfread(chanFile);
%channels = struct2table(channels);

% Read in the events file
eventsFile = fullfile(sessionDir, 'ieeg', sprintf('%s_events.tsv', fname_in));
if ~exist(eventsFile, 'file'), error('events file not found: %s', eventsFile); end
%fprintf('[%s] Reading in events file: %s\n', mfilename, eventsFile); 
fprintf('.'); 
events = readtable(eventsFile, 'FileType', 'text');
%events = tdfread(eventsFile);
%events = struct2table(events);
 
% Read in the data file
dataFile = fullfile(sessionDir, 'ieeg', sprintf('%s_ieeg.eeg', fname_in));
if ~exist(dataFile, 'file'), error('data file not found: %s', dataFile); end
%fprintf('[%s] Reading in data file: %s\n', mfilename, dataFile); 
fprintf('.'); 
hdr = ft_read_header(dataFile);
data = ft_read_data(dataFile);

% Read in the json file
jsonReadFile = fullfile(sessionDir, 'ieeg', sprintf('%s_ieeg.json', fname_in));
if ~exist(jsonReadFile, 'file'), error('data file not found: %s', dataReadFile); end
%fprintf('[%s] Reading in ieeg json file: %s\n', mfilename, jsonReadFile); 
fprintf('.\n'); 
ieeg_json = jsonread(jsonReadFile);
           
                           
end

