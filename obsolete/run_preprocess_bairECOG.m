% run_preprocess_bairECOG

% SCRIPT DESCRIPTION %

% This script calls the preprocessing script for the BIDS-formatted BAIR
% ECoG data for a specific subject and session. Data should be in BIDS.

%% Define dataset to preprocess %%

% Dataset specs
project_name = 'visual';

sub_label    = 'som748';
ses_label   = {'nyuecog02'};

%sub_label   = 'som726'; 
%ses_label   = {'nyuecog01','nyuecog02','nyuecog03'};

%sub_label   = 'umcuchaam'; 
%ses_label   = {'UMCUECOGday03'};

%sub_label   = 'beilen'; 
%ses_label   = {'UMCUECOGday03', 'UMCUECOGday06'};

%sub_label   = 'intraop015'; 
%ses_label   = {'UMCUOR'};

% PREPROCESSING SPECS

% broadband
%specs.bb_bands     = [[70 90]; [90 110]; [130 150]; [150 170]];
%specs.bb_bands     = [[60 70]; [70 80]; [80 90]; [110 120]; [120 130]; [130 140]; [160 170]];    
specs.bb_bands     = [[70 80]; [80 90]; [90 100]; [100 110]; [130 140]; [140 150]; [150 160]; [160 170]];
%specs.bb_bands     = [[70 80]; [80 90]; [90 100]; [100 110]; [130 140]; [140 150]; [150 160]; [160 170]; [190 200]];
specs.bb_method    = [];
specs.bb_bands     = [[70 80]; [80 90]; [90 100]; [100 110]; [130 140]];

% epochs
specs.epoch        = [-0.5 1.5]; % seconds
% specs.epoch      = [-0.5 1.5]; % Kenichi

% general
specs.make_plots   = 'yes'; % 'yes', 'no'
specs.overwrite    = 'yes'; % 'yes', 'no'

%% Run!

dataPth = fullfile('/Volumes/server/Projects/BAIR/Data/BIDS/', project_name);

% Check whether we have the ECoG_utils repository on the path
if ~exist('createBIDS_ieeg_json_nyuSOM.m')
    tbUse ECoG_utils;
end

for ii = 1:length(ses_label)
    preprocess_bairECOG(dataPth, sub_label, ses_label{ii}, specs);
end




