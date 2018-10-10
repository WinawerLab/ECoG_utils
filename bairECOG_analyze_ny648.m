tbUse('ECoG_utils');

% SCRIPT DESCRIPTION %

%% [0] Define paths and dataset specs %%

% Dataset specs
projectName = 'visual';
sub_label   = 'som648'; 
ses_label   = 'nyuecog01';

% Input paths specs
dataPth     = '/Volumes/server/Projects/BAIR/Data/BIDS/';

% Output paths specs
dataDir = fullfile(dataPth, projectName, 'derivatives', 'preprocessed', sprintf('sub-%s', sub_label));

%% Load preprocessed data

dataName = fullfile(dataDir, sprintf('sub-%s_ses-%s_preproc', sub_label, ses_label));
load(dataName);

dataName = fullfile(dataDir, sprintf('sub-%s_ses-%s_epoched', sub_label, ses_label));
load(dataName);

%% Plot responses for visual electrodes

% plot specs
whichElectrodes = trials.viselec.benson14_varea.elec_labels;
trialType = {'CRF','ONEPULSE', 'TWOPULSE'}; %{'HRFPATTERN'};
%trialType = {'HRFPATTERN'};
%trialType = {'SCENES', 'FACES', 'LETTERS'};
collapseTrialTypes = 'no';

% Plot both broadband and evoked response per electrode
ecog_plotTrials(trials, whichElectrodes, trialType, collapseTrialTypes);

