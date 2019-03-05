% run_preprocess_bairECOG

% SCRIPT DESCRIPTION %

% This script calls the preprocessing script for the BIDS-formatted BAIR
% ECoG data for a specific subject and session.

% [0] Path specifications: data should be in BIDS

%% Define dataset to preprocess %%

% Dataset specs

project_name = 'visual';

% sub_label   = 'chaam'; 
% ses_label   = 'UMCUECOGday03';

sub_label   = 'som708'; 
ses_label   = {'nyuecog01','nyuecog02'};

% PREPROCESSING SPECS

% broadband
%specs.bb_bands     = [[70 90]; [90 110]; [130 150]; [150 170]];
%specs.bb_bands     = [[60 70]; [70 80]; [80 90]; [110 120]; [120 130]; [130 140]; [160 170]];    
specs.bb_bands     = [[70 80]; [80 90]; [90 100]; [100 110]; [130 140]; [140 150]; [150 160]; [160 170]];
specs.bb_method    = 7;

% epochs
specs.epoch        = [-0.5 1.5]; % seconds

% baseline correction

% spectra
specs.reg_erp      = 0; % regress out erp for spectra?

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




