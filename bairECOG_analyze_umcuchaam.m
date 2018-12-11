
% SCRIPT DESCRIPTION %

%% [0] Define paths and dataset specs %%

% Dataset specs
projectName = 'visual';
sub_label   = 'chaam'; 
ses_label   = 'UMCUECOGday03';

% Input paths specs
dataPth     = '/Volumes/server/Projects/BAIR/Data/BIDS/';

% Output paths specs
procDir = fullfile(dataPth, projectName, 'derivatives', 'preprocessed', sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label));

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
dataName = fullfile(procDir, sprintf('sub-%s_ses-%s_epoched.mat', sub_label, ses_label));
load(dataName);

%% Do electrode matching
dataDir = fullfile(dataPth, projectName, sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label), 'ieeg');

specs = [];
specs.pID           = sub_label; % patient ID number
specs.patientPool   = 'BAIR';
specs.fsDir         = '/Volumes/server/Freesurfer_subjects/umcuchaam';
specs.plotmesh      = 'right';
specs.plotelecs     = 'yes';
trials.viselec      = electrode_to_nearest_node(specs, dataDir);

%% Plot responses for visual electrodes

% FROM GIO/DORA: Electrodes on seizure are OC12, OC13, OC21, OC22.
ElecsToExclude  = [trials.channels.name(contains(trials.channels.status, 'bad'))' {'Oc9', 'Oc12', 'Oc13', 'Oc21', 'Oc22'}]; % Oc9 labeled as bad

% Which electrodes to plot? (Each electrode gets a subplot)
%whichElectrodes = unique([trials.viselec.benson14_varea.elec_labels trials.viselec.wang15_mplbl.elec_labels]); 
%whichElectrodes = setdiff(whichElectrodes,ElecsToExclude);

% SELECTING ONLY V1, V2, V3, V3a, V3b:
whichElectrodes = trials.viselec.benson14_varea.elec_labels(strncmp(trials.viselec.benson14_varea.area_labels, 'V',1));
whichElectrodes = ecog_sortElectrodesonVisualArea(whichElectrodes,trials.viselec.benson14_varea);

% Which stimulus conditions to plot? 
%whichTrials = {}; % if empty, plot average of all trials
%whichTrials = {'CRF','SPARSITY','ONEPULSE', 'TWOPULSE', 'GRATING', 'PLAID', 'CIRCULAR', 'HOUSES', 'FACES', 'LETTERS'}; % everything except prf
%collapseTrialTypes = 'yes'; % if yes, plots average across all trial types listed under 'whichTrials'

%whichTrials = {'CRF-1','CRF-2', 'CRF-3','CRF-4', 'CRF-5', 'HRFPATTERN'};
%whichTrials = {'ONEPULSE-1','ONEPULSE-2', 'ONEPULSE-3','ONEPULSE-4', 'ONEPULSE-5', 'ONEPULSE-6'};
whichTrials = {'ONEPULSE-5', 'TWOPULSE-1','TWOPULSE-2', 'TWOPULSE-3','TWOPULSE-4', 'TWOPULSE-5', 'TWOPULSE-6'};

%whichTrials = {'HOUSES','FACES', 'LETTERS'};
%whichTrials = {'GRATING', 'PLAID','CIRCULAR'};
%whichTrials = {'CRF', 'SPARSITY','ONEPULSE','GRATING', 'PLAID','CIRCULAR'};
%whichTrials = {'SPARSITY-1','SPARSITY-2', 'SPARSITY-3','SPARSITY-4'};
%whichTrials = {'CRF-5','ONEPULSE-5', 'TWOPULSE-5'}; 
%whichTrials = {'HRF'};
%collapseTrialTypes = 'no'; 

%baselineType = 'selectedtrials';%'selectedtrials';

% Smooth the timecourses?
% defined in milliseconds; if 0 or empty, no smoothing will be applied
%smoothingLevelInMs = [20]; 

% Plot both broadband and evoked response per electrode
% 'out' contains the plotted time courses and SEs
[out] = ecog_plotTrials(trials, whichElectrodes, whichTrials);



