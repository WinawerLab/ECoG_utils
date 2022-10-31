function [vertices_r, faces_r, vertices_l, faces_l]= bidsEcogReadSurfFile(dataPath, subject, surfaceType)

% Reads in freesurfer surface coordinate file from dataPath with BIDS formatted data. 
%
% Input
%     dataPath:         path where the BIDS project lies (string)
%     subject:          BIDS subject name (string, all lower case)
%     surfaceType:      freesurfer surface type (string). default 'pial'
%
% Output
%     [vertices_r, faces_r, vertices_l, faces_l]
%
% Example
%   dataPath = '/Volumes/server/Projects/BAIR/Data/BIDS/visual';
%   subject = 'p10';
%  [vertices_r, faces_r, vertices_l, faces_l]= bidsEcogReadSurfFile(projectDir, subject)
%
% Makes use of 'read_surf' function from freesurfer
%
% IG, BAIR 2020

% <dataPath>
if ~exist('dataPath', 'var') || isempty(dataPath), error('dataPath not defined'); end    

% <subject>
if ~exist('subject', 'var') || isempty(subject), error('subject not defined'); end

% <surfaceType>
if ~exist('surfaceType', 'var') || isempty(surfaceType), surfaceType = 'pial'; end

% Define path
fsDir = fullfile(dataPath, 'derivatives', 'freesurfer',  sprintf('sub-%s', subject));
%fprintf('[%s] Reading from %s\n', mfilename, fsDir);

% Define filenames
surf_file_rh = fullfile(fsDir, 'surf', sprintf('rh.%s', surfaceType));
surf_file_lh = fullfile(fsDir, 'surf', sprintf('lh.%s', surfaceType));

% Read in the files
if exist(surf_file_rh, 'file') && exist(surf_file_lh, 'file')
    fprintf('[%s] Reading Freesurfer %s surface file %s ...\n',mfilename, surfaceType, surf_file_rh);
    [vertices_r, faces_r] = read_surf(surf_file_rh);
    fprintf('[%s] Reading Freesurfer %s surface file %s ...\n',mfilename, surfaceType, surf_file_lh);
    [vertices_l, faces_l] = read_surf(surf_file_lh);
else
    error('No freesurfer surfaces found for sub-%s in %s',subject, fsDir);
end

end