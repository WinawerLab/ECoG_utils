function [electrode_table] = bidsEcogReadElectrodeFile(dataPath, subject, session)

% Reads in ECoG electrode coordinate file from dataPath with BIDS formatted data. 
%
% Input
%     dataPath:         path where the BIDS project lies (string)
%     subject:          BIDS subject name (string, all lower case)
%     session:          BIDS session name (string, all lower case)
%
% Output
%     electrode_table:  A Matlab table with the contents of electrodes.tsv
%
% Example
%   dataPath = '/Volumes/server/Projects/BAIR/Data/BIDS/visual';
%   subject = 'som726';
%   session = 'nyuecog01'
%  [electrode_table, elec_xyz] = bidsEcogReadElectrodeFile(dataPath, subject, session)
%
% See also bidsEcogReadFiles.m bidsSpecifySessions.m
%
% IG, BAIR 2020


% <dataPath>
if ~exist('dataPath', 'var') || isempty(dataPath), error('dataPath not defined'); end     

% <subject>
if ~exist('subject', 'var') || isempty(subject), error('subject not defined'); end

% <session>
if ~exist('session', 'var') || isempty(session), error('session not defined'); end

sessionDir = fullfile(dataPath, sprintf('sub-%s', subject), sprintf('ses-%s', session), 'ieeg');
fprintf('[%s] Reading from %s\n', mfilename, sessionDir);

D = dir(fullfile(sessionDir, '*electrodes.tsv'));          
if isempty(D)
    error('No electrode coordinate file found in %s', sessionDir);
	%fprintf('[%s] Electrode coordinate file not found in %s - exiting [%s]. \n', sessionDir);
    %electrode_table = [];
else
    elec_file = D(1).name;
    fprintf('[%s] Reading %s...\n',mfilename,fullfile(sessionDir,elec_file));
    electrode_table = readtable(fullfile(sessionDir, elec_file), 'FileType', 'text');
end

end
