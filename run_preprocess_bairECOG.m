% run_preprocess_bairECOG

% SCRIPT DESCRIPTION %

% This script calls the preprocessing script for the BIDS-formatted BAIR
% ECoG data for a specific subject and session.

% [0] Path specifications: data should be in BIDS

%% Define dataset to preprocess %%

% Dataset specs

project_name = 'visual';
sub_label   = 'som648'; 
ses_label   = 'nyuecog01';

% Preprocessing specs

specs.bb_bands     = [[70 90]; [90 110]; [130 150]; [150 170]];
specs.bb_method    = 7;
specs.epoch        = [-0.5 1.5]; % seconds
specs.make_plots   = 'yes'; % 'yes', 'no'
specs.overwrite    = 'no'; % 'yes', 'no'

%% Run!

dataPth = fullfile('/Volumes/server/Projects/BAIR/Data/BIDS/', project_name);

% Check whether we have the ECoG_utils repository on the path
if ~exist('createBIDS_ieeg_json_nyuSOM.m')
    tbUse ECoG_utils;
end

preprocess_bairECOG(dataPth, sub_label, ses_label, specs);