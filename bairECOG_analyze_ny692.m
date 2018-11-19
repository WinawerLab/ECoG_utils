
% SCRIPT DESCRIPTION %

%% [0] Define paths and dataset specs %%

% Dataset specs
projectName = 'visual';
sub_label   = 'som692'; 
ses_label   = 'nyuecog01';

% Input paths specs
dataPth     = '/Volumes/server/Projects/BAIR/Data/BIDS/';

% Output paths specs
dataDir = fullfile(dataPth, projectName, 'derivatives', 'preprocessed', sprintf('sub-%s', sub_label));

%% Load preprocessed data

% This file contains the full concatenated data (raw voltages, common
% average rereferenced voltages, broadband time series):
%dataName = fullfile(dataDir, sprintf('sub-%s_ses-%s_preproc', sub_label, ses_label));
%load(dataName);

% This file contains stimulus-defined epochs of the rereferenced voltages
% (trials.evoked) and the broadband time series (trials.broadband) plus a
% separate trials struct which has 'baseline' epochs (no stimulus shown).
% NOTE that prestimulus baseline correction has not yet been performed for
% the time courses, this now happens in the plotting script below (should
% probably be a separate step)
dataName = fullfile(dataDir, sprintf('sub-%s_ses-%s_epoched', sub_label, ses_label));
load(dataName);

%% Plot responses for visual electrodes

% Which electrodes to plot? (Each electrode gets a subplot)
% e.g. 'MO01'; e.g. trials.viselec.benson14_varea.elec_labels
%whichElectrodes = unique([trials.viselec.benson14_varea.elec_labels trials.viselec.wang2015_atlas.elec_labels]); 
%whichElectrodes = {'MO01', 'MO02','MO03','MO04'};
whichElectrodes = {'BO01','BO02','BO03','BO04','BT01','BT02','BT03','BT04','G03','G06','G07','G08','G16','G23','G24','G31','G32','MIO01','MIO02','MIO03','MIO04','MSO02','MSO03','MSO04'};
%whichElectrodes = {'BO01','BO02','BO03','BO04','BT01','G07','G08','G24','MIO01','MIO02','MIO03','MIO04'};
%whichElectrodes = {'BO01','BO02','MIO02','BO04'};
%whichElectrodes = {'BO01'};

% Which stimulus conditions to plot? 
%whichTrials = {}; % if empty, plot average of all trials
%whichTrials = {'CRF','SPARSITY','ONEPULSE', 'TWOPULSE', 'GRATING', 'PLAID', 'CIRCULAR', 'HOUSES', 'FACES', 'LETTERS'}; % everything except prf
%collapseTrialTypes = 'yes'; % if yes, plots average across all trial types listed under 'whichTrials'

%whichTrials = {'HOUSES','FACES', 'LETTERS'};
%whichTrials = {'GRATING', 'PLAID','CIRCULAR'};
%whichTrials = {'CRF', 'SPARSITY','ONEPULSE','GRATING', 'PLAID','CIRCULAR'};
%whichTrials = {'CRF-1','CRF-2', 'CRF-3','CRF-4', 'CRF-5'};
whichTrials = {'SPARSITY-1','SPARSITY-2', 'SPARSITY-3','SPARSITY-4'};
%whichTrials = {'ONEPULSE-1','ONEPULSE-2', 'ONEPULSE-3','ONEPULSE-4', 'ONEPULSE-5'};
%whichTrials = {'TWOPULSE-1','TWOPULSE-2', 'TWOPULSE-3','TWOPULSE-4', 'TWOPULSE-5'};
%whichTrials = {'CRF-5','ONEPULSE-5', 'TWOPULSE-5'}; 
collapseTrialTypes = 'no'; 

baselineType = 'selectedtrials';%'selectedtrials';

% Smooth the timecourses?
% defined in milliseconds; if 0 or empty, no smoothing will be applied
smoothingLevelInMs = [20]; 

% Plot both broadband and evoked response per electrode
% 'out' contains the plotted time courses and SEs
[out] = ecog_plotTrials(trials, whichElectrodes, whichTrials, collapseTrialTypes, smoothingLevelInMs, baselineType);



