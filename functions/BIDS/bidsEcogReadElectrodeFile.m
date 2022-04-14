function [electrode_table, elec_xyz] = bidsEcogReadElectrodeFile(dataPath, subject, session)

% Reads in ECoG electrode coordinate file from dataPath with BIDS formatted data. 
%
% Input
%     dataPath:         path where the BIDS project lies (string)
%     subject:          BIDS subject name (string, all lower case)
%     session:          BIDS session name (string, all lower case)
%
% Output
%     electrode_table:  A Matlab table with the contents of electrodes.tsv
%     elec_xyz:         Electrode locations with shift based on 
%                       electrodes_shift.tsv in freesurfer directory
%
% Example
%   dataPath = '/Volumes/server/Projects/BAIR/Data/BIDS/visual_ecog_recoded';
%   subject = 'p10';
%   session = 'nyuecog01'
%  [electrode_table, elec_xyz] = bidsEcogReadElectrodeFile(dataPath, subject, session)
%
% See also bidsEcogReadFiles.m bidsSpecifySessions.m
%
% IG, BAIR 2020; K.Yuasa, 2022


% <dataPath>
if ~exist('dataPath', 'var') || isempty(dataPath), error('dataPath not defined'); end     

% <subject>
if ~exist('subject', 'var') || isempty(subject), error('subject not defined'); end

% <session>
if ~exist('session', 'var') || isempty(session), error('session not defined'); end

sessionDir = fullfile(dataPath, sprintf('sub-%s', subject), sprintf('ses-%s', session), 'ieeg');
fsDir      = fullfile(dataPath, 'derivatives', 'freesurfer',  sprintf('sub-%s', subject));
%fprintf('[%s] Reading from %s\n', mfilename, sessionDir);

% Load electrode coordinate file
D = dir(fullfile(sessionDir, '*electrodes.tsv'));          
if isempty(D)
    error('No electrode coordinate file found in %s', sessionDir);
	%fprintf('[%s] Electrode coordinate file not found in %s - exiting [%s]. \n', sessionDir);
    %electrode_table = [];
else
    elec_file = D(1).name;
    fprintf('[%s] Reading electrode file %s...\n',mfilename,fullfile(sessionDir,elec_file));
    electrode_table = readtable(fullfile(sessionDir, elec_file), 'FileType', 'text');
end

% Set electrode locations
elec_xyz    = [electrode_table.x electrode_table.y electrode_table.z];

% In some patients (e.g. UMCU patients) the center coordinate of the brain
% image doesn't match to the electrode coordinates. The amount of shift
% needed is provided in a custom-made file 'electrode_shift.tsv' that
% we place inside the freesurfer mri directory. The file has 'c_r', 'c_a',
% and 'c_s' in the first row, and the shift amount for each in the next. 
% The values in there are based on the ras xform info inside the orig.mgz
% header in the freesurfer mri directory (to be inspected using
% freesurfer's mri_info orig.mgz)

% Shift electrode locations
elec_shift_file = fullfile(fsDir, 'mri', 'electrode_shift.tsv');
if exist(elec_shift_file,'file')
    elec_shift = readtable(elec_shift_file, 'FileType', 'text');
    if istablefield(elec_shift,'c_r')
        elec_xyz(:,1) = elec_xyz(:,1) - elec_shift.c_r;
    end
    if istablefield(elec_shift,'c_a')
        elec_xyz(:,2) = elec_xyz(:,2) - elec_shift.c_a;
    end
    if istablefield(elec_shift,'c_s')
        elec_xyz(:,3) = elec_xyz(:,3) - elec_shift.c_s;
    end    
end
