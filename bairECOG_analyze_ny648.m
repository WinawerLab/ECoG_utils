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
procDir = fullfile(dataPth, projectName, 'derivatives', 'preprocessed', sprintf('sub-%s', sub_label));

%% Load preprocessed data

% This file contains the full concatenated data (raw voltages, common
% average rereferenced voltages, broadband time series):
% dataName = fullfile(dataDir, sprintf('sub-%s_ses-%s_preproc', sub_label, ses_label));
% load(dataName);

% This file contains stimulus-defined epochs of the rereferenced voltages
% (trials.evoked) and the broadband time series (trials.broadband) plus a
% separate trials struct which has 'baseline' epochs (no stimulus shown).
% NOTE that prestimulus baseline correction has not yet been performed for
% the time courses, this now happens in the plotting script below (should
% probably be a separate step)
dataName = fullfile(procDir, sprintf('sub-%s_ses-%s_epoched', sub_label, ses_label));
load(dataName);


%% compute spectra
% for ii = size(trials.channels,1) 
%     [Pxx,f] = pwelch(squeeze(trials.evoked(ii,:,:)),trials.fsample,0,trials.fsample,trials.fsample);
%     trials.spectra(ii,:,:) = Pxx;  
% end
% trials.f = f;

%% plot electrodes on mesh
dataDir = fullfile(dataPth, projectName, sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label), 'ieeg');

specs = [];
specs.pID           = sub_label; % patient ID number
specs.patientPool   = 'BAIR';
specs.plotmesh      = 'both';
visualelectrodes    = electrode_to_nearest_node(specs, dataDir);

%% Plot responses for visual electrodes

% Which electrodes to plot? (Each electrode gets a subplot)
% e.g. 'MO01'; e.g. trials.viselec.benson14_varea.elec_labels
%whichElectrodes = unique([trials.viselec.benson14_varea.elec_labels, trials.viselec.wang2015_atlas.elec_labels]); 
%whichElectrodes = {'IO01', 'IO02','IO03','IO04'};
whichElectrodes = {'MO01', 'MO02','MO03','MO04'};
%whichElectrodes = {'IO02'};
%whichElectrodes = {'MO01'};

% Which stimulus conditions to plot? 
%whichTrials = {'HRFPATTERN','CRF','ONEPULSE', 'TWOPULSE'}; 
%whichTrials = {'HRFPATTERN', 'SCENES','FACES', 'LETTERS'};
%whichTrials = {'CIRCULAR', 'GRATING','PLAID'};
%whichTrials = {'HRFPATTERN','CRF-1','CRF-2', 'CRF-3','CRF-4', 'CRF-5'};
%whichTrials = {'SPARSITY-1','SPARSITY-2', 'SPARSITY-3','SPARSITY-4'};
%whichTrials = {'ONEPULSE-1','ONEPULSE-2', 'ONEPULSE-3','ONEPULSE-4', 'ONEPULSE-5'};
whichTrials = {'TWOPULSE-1','TWOPULSE-2', 'TWOPULSE-3','TWOPULSE-4', 'TWOPULSE-5'};
%whichTrials = {'HRFPATTERN','CRF-5','ONEPULSE-5', 'TWOPULSE-5'}; 
%whichTrials = {}; % if empty, plot average of all trials
%whichTrials = {'HRFPATTERN'};

% Plot each type individually or combine?
% if yes, plots average across all trial types listed under 'whichTrials'
collapseTrialTypes = 'no'; 

% Smooth the timecourses?
% defined in milliseconds; if 0 or empty, no smoothing will be applied
smoothingLevelInMs = []; 

% Plot both broadband and evoked response per electrode
% 'out' contains the plotted time courses and SEs
[out] = ecog_plotTrials(trials, whichElectrodes, whichTrials, collapseTrialTypes, smoothingLevelInMs);



