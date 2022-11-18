function [curvs_r, curvs_l]= bidsEcogReadCurvFile(dataPath, subject, curvType)

% Reads in freesurfer surface coordinate file from dataPath with BIDS formatted data. 
%
% Input
%     dataPath:         path where the BIDS project lies (string)
%     subject:          BIDS subject name (string, all lower case)
%     curvType:         freesurfer curvature type (string). default 'sulc'
%
% Output
%     [curvs_r, curvs_l]
%
% Example
%   dataPath = '/Volumes/server/Projects/BAIR/Data/BIDS/visual';
%   subject = 'p10';
%  [curvs_r, curvs_l]= bidsEcogReadCurvFile(projectDir, subject)
%
% Makes use of 'read_curv' function from freesurfer
%
% KY, BAIR 2022

% <dataPath>
if ~exist('dataPath', 'var') || isempty(dataPath), error('dataPath not defined'); end    

% <subject>
if ~exist('subject', 'var') || isempty(subject), error('subject not defined'); end

% <surfaceType>
if ~exist('curvType', 'var') || isempty(curvType), surfaceType = 'sulc'; end

% Define path
fsDir = fullfile(dataPath, 'derivatives', 'freesurfer',  sprintf('sub-%s', subject));
%fprintf('[%s] Reading from %s\n', mfilename, fsDir);

% Define filenames
curv_file_rh = fullfile(fsDir, 'surf', sprintf('rh.%s', curvType));
curv_file_lh = fullfile(fsDir, 'surf', sprintf('lh.%s', curvType));

% Read in the files
if exist(curv_file_rh, 'file') && exist(curv_file_lh, 'file')
    fprintf('[%s] Reading Freesurfer %s curvature file %s ...\n',mfilename, curvType, curv_file_rh);
    [curvs_r] = read_curv(curv_file_rh);
    fprintf('[%s] Reading Freesurfer %s curvature file %s ...\n',mfilename, curvType, curv_file_lh);
    [curvs_l] = read_curv(curv_file_lh);
else
    error('No freesurfer curvatures found for sub-%s in %s',subject, fsDir);
end

end