function [coordsystem_json, json_options] = createBIDS_ieeg_coordsystem_json_nyuSOM()

% Adapted from BIDS_starter_kit createBIDS_ieeg_coordystem_json template
% Intended to be run to get default options for NYU SOM recordings
%
% IGroen, 2018

%% Template Matlab script to create an BIDS compatible _coordsystem.json file
% For BIDS-iEEG
% This example lists all required and optional fields.
% When adding additional metadata please use CamelCase 
%
% DHermes, 2017
% modified RG 201809

%%  Required fields

coordsystem_json.iEEGCoordinateSystem  = 'Other'; % REQUIRED. Defines the coordinate 
% system for the iEEG electrodes. See Appendix VIII: preferred names of Coordinate 
% systems. If "Other" (e.g. individual subject MRI), provide definition of the 
% coordinate system in  [iEEGCoordinateSystemDescription]. If positions correspond to 
% pixel indices in a 2D image (of either a volume-rendering, surface-rendering, 
% operative photo, or operative drawing), this must be ?pixels?. 

coordsystem_json.iEEGCoordinateUnits  = 'mm'; %REQUIRED. Units of the _electrodes.tsv, 
% MUST be ?m?, ?mm?, ?cm? or ?pixels?. 

%% Optional fields

coordsystem_json.iEEGCoordinateSystemDescription  = 'T1';% %OPTIONAL. Freeform text description or 
% link to document describing the iEEG coordinate system system in detail (e.g. ?ACPC?). 

coordsystem_json.IntendedFor = ''; % REQUIRED. This can be an MRI/CT or a file containing 
% the operative photo, x-ray or drawing with path relative to the project folder. If only 
% a surface reconstruction is available, this should point to the surface reconstruction file. 
% Note that this file should have the same coordinate system specified in iEEGCoordinateSystem. 

coordsystem_json.iEEGCoordinateProcessingDescription = 'surface_projection'; % REQUIRED. Has any 
% projection been done on the electrode positions (e.g. ?surface_projection?,  ?none?).

coordsystem_json.iEEGCoordinateProcessingReference = 'PMID: 22759995'; % RECOMMENDED. A reference to a paper that defines in more detail 
% the method used to project or localize the electrodes

%% 
json_options.indent = '    '; % this just makes the json file look prettier 
% when opened in a text editor

end
