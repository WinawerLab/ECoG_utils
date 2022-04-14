function [atlas_r, atlas_l]= bidsEcogReadAtlasFile(dataPath, subject, atlasName)

% Reads in freesurfer surface coordinate file from dataPath with BIDS formatted data. 
%
% Input
%     dataPath:         path where the BIDS project lies (string)
%     subject:          BIDS subject name (string, all lower case)
%     atlasName:        atlas filename (string), e.g. 'wang15_mplbl'
%
% Output
%     [atlas_r, atlas_l]
%
% Example
%   dataPath  = '/Volumes/server/Projects/BAIR/Data/BIDS/visual';
%   subject   = 'p10';
%   atlasName = 'benson14_varea';
%  [atlas_r, atlas_l]= bidsEcogReadAtlasFile(projectDir, subject, atlasName)
%
% Makes use of 'load_mgh' function from vistasoft
%
% IG, BAIR 2020

% <dataPath>
if ~exist('dataPath', 'var') || isempty(dataPath), error('dataPath not defined'); end    

% <subject>
if ~exist('subject', 'var') || isempty(subject), error('subject not defined'); end

% <surfaceType>
if ~exist('atlasName', 'var') || isempty(atlasName), error('atlasName not defined'); end

% Define path
fsDir = fullfile(dataPath, 'derivatives', 'freesurfer',  sprintf('sub-%s', subject));
fprintf('[%s] Reading atlas file %s from %s\n', mfilename, atlasName, fsDir);

% Define filenames
atlas_file_rh = fullfile(fsDir, 'surf', sprintf('rh.%s.mgz', atlasName));
atlas_file_lh = fullfile(fsDir, 'surf', sprintf('lh.%s.mgz', atlasName));

% Read in the files
if exist(atlas_file_rh, 'file') && exist(atlas_file_lh, 'file')
    [atlas_r] = load_mgh(atlas_file_rh);
    [atlas_l] = load_mgh(atlas_file_lh);
else
    error('Atlas %s not found',atlasName);
    %atlas_r = [];
    %atlas_l = [];
end

end