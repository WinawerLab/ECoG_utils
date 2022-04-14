function [fname_out] = bidsEcogWriteFiles(dataPath, subject, session, task, runnum, description, ...
    data, channels, events, ieeg_json, hdr)
% Writes out ECoG data files to dataPath in BIDS format.
%
% [fname_out] = bidsEcogWriteFiles(dataPath, subject, session, task, ...
%    runnum, description, data, channels, events, ieeg_json, hdr)
%
% Input
%     dataPath:         path where the BIDS project lies (string)
%     subject:          BIDS subject name (string, all lower case)
%     session:          BIDS session name (string, all lower case)
%     task:             BIDS task name (string, all lower case)                           
%     runnum:           BIDS run number (string, all lower case)
%     description:      String for the BIDS 'desc-' label in the bidsname.
%                       To omit, leave empty ([]).
%     data:             A 2D matrix with ECoG data (channels x time)
%     channels:         A table with channel information
%     events:           A table with events information
%     ieeg_json:        A struct with the ECoG meta data
%     hdr:              The data header as provided by Fieldtrip (required
%                       in order to write data files)
%
% Output
%     fname_out:        The bidsname used to name the data files.
%
% Note
% This function uses FieldTrips ft_write_data.m to write out the data,
% which is why the Fieldtrip header needs to provided as input. The data
% will be written out in BVA (brain vision analyzer) format, i.e. as .eeg,
% .vmrk and .vhdr. In addition the following meta-data files will be
% written: a channel file (<bidsfilename>_channels.tsv), an events file
% (<bidsfilename>_events.tsv) and a json file (<bidsfilename>_ieeg.json). 
%
% Example
%
% % Read in data:
% projectDir = '/Volumes/server/Projects/BAIR/Data/BIDS/visual';
% dataPath = fullfile(projectDir, 'derivatives', 'ECoGBroadband');
% subject = 'p10';
% session = 'nyuecog01'
% task = 'prf';
% runnum = '01';
% description = 'broadband';
% [data, channels, events, ieeg_json, hdr] = bidsEcogReadFiles(dataPath, ...
% subject, session, task, runnum, description);
%
% % Write out data:
% [fname_out] = bidsEcogWriteFiles(dataPath, subject, session, task, ...
%    runnum, description, data, channels, events, ieeg_json, hdr);
% 
% See also bidsEcogReadFiles.m bidsSpecifyData.m bidsSpecifySessions.m
% See also ft_write_data.m 

% <description>
if ~exist('description', 'var')
    description = [];
end

%fprintf('[%s] Task = %s, Run = %s \n', mfilename, task, runnum);

if ~isempty(description)
    fname_out = sprintf('sub-%s_ses-%s_task-%s_run-%s_desc-%s', subject, session, task, runnum, description);
else
    fname_out = sprintf('sub-%s_ses-%s_task-%s_run-%s', subject, session, task, runnum);
end

sessionDir = fullfile(dataPath, sprintf('sub-%s', subject), sprintf('ses-%s', session), 'ieeg');
if ~exist(sessionDir, 'dir'), mkdir(sessionDir);end
fprintf('[%s] Writing to %s\n', mfilename, sessionDir);
fprintf('[%s] Writing %s', mfilename, fname_out);

% Write out the channels file
chanWriteFile = fullfile(sessionDir, sprintf('%s_channels.tsv', fname_out));
%fprintf('[%s] Writing new channels file: %s\n', mfilename, chanWriteFile); 
fprintf('.'); 
writetable(channels,chanWriteFile,'FileType','text','Delimiter','\t');
       
% Write out the events file
eventsWriteFile = fullfile(sessionDir, sprintf('%s_events.tsv', fname_out));
%fprintf('[%s] Writing new channels file: %s\n', mfilename, eventsWriteFile); 
fprintf('.'); 
writetable(events,eventsWriteFile,'FileType','text','Delimiter','\t');
  
% Save out the data file
dataWriteFile = fullfile(sessionDir, sprintf('%s_ieeg.eeg', fname_out));
%fprintf('[%s] Writing new data file: %s\n', mfilename, dataWriteFile); 
fprintf('.'); 
ft_write_data(dataWriteFile, data, 'header', hdr, 'dataformat', 'brainvision_eeg');   

% Save out the ieeg_json file
jsonWriteFile = fullfile(sessionDir, sprintf('%s_ieeg.json', fname_out));
%fprintf('[%s] Writing new ieeg json file: %s\n', mfilename, jsonWriteFile); 
fprintf('.\n'); 
json_options.indent = '    '; 
jsonwrite(jsonWriteFile,ieeg_json, json_options)
                                      
end

