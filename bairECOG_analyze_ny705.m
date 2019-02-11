
% SCRIPT DESCRIPTION %

%% [0] Define paths and dataset specs %%

% Dataset specs
projectName = 'motor';
sub_label   = 'ny705'; 
ses_label   = 'nyuecog01';

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
dataName = fullfile(procDir, sprintf('sub-%s_ses-%s_epoched', sub_label, ses_label));
load(dataName);

% %% plot electrodes on mesh
% dataDir = fullfile(dataPth, projectName, sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label), 'ieeg');
% 
% specs = [];
% specs.pID           = sub_label; % patient ID number
% specs.patientPool   = 'BAIR';
% specs.plotmesh      = 'both';
% visualelectrodes    = electrode_to_nearest_node(specs, dataDir);

%% Plot responses 

% Which electrodes to plot? (Each electrode gets a subplot)
% e.g. 'MO01'; e.g. trials.viselec.benson14_varea.elec_labels
whichElectrodes = trials.channels.name(contains(trials.channels.name, 'GB'));
%whichElectrodes = {'GB23','GB35','GB36','GB37', 'GB38' 'GB53', 'GB69','GB70','GB71','GB85','GB86','GB87'};

% Which stimulus conditions to plot? 
%whichTrials = {}; % if empty, plot average of all trials
%whichTrials = {'CRF','SPARSITY','ONEPULSE', 'TWOPULSE', 'GRATING', 'PLAID', 'CIRCULAR', 'HOUSES', 'FACES', 'LETTERS'}; % everything except prf

whichTrials = {'CLENCH'};
%whichTrials = {'D', 'Y', 'V', 'F'};
%whichTrials = {'thumb','index','middle','ring','pinky'};

% Plot both broadband and evoked response per electrode
% 'out' contains the plotted time courses and SEs

% PLOT time course for each condition
specs = [];
specs.dataTypes          = {'broadband'};
specs.smoothingLevelInMs = [];
specs.collapseTrialTypes = 'yes';
specs.baselineType       = 'selectedtrials';%'selectedtrials';
specs.plot.colorMap      = 'parula';
specs.plot.nSubPlots     = [];
specs.plot.addEccToTitle = 'yes';
specs.plot.showMax       = 'no';

[out] = ecog_plotTimecourses(trials, whichElectrodes, whichTrials, specs);

%elecInx = 1:63:length(whichElectrodes)+63;
% for ii = 1:length(elecInx)-1
%     [out] = ecog_plotTimecourses(trials, whichElectrodes(elecInx(ii):elecInx(ii+1)-1), whichTrials, specs);
% end



