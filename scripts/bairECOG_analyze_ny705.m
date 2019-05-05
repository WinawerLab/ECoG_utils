
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

%% plot electrodes on mesh
dataDir = fullfile(dataPth, projectName, sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label), 'ieeg');

specs = [];
specs.pID           = 'som705'; % patient ID number
specs.atlasNames    = {'benson14_varea'};
specs.patientPool   = 'BAIR';
specs.plotmesh      = 'both';
specs.plotlabel     = 'no';
visualelectrodes    = electrode_to_nearest_node(specs, dataDir);

%% Plot responses 

% Which electrodes to plot? (Each electrode gets a subplot)
whichElectrodes = trials.channels.name(contains(trials.channels.name, 'GB'));
whichElectrodes = [{'GB1-nodata'};{'GB2-nodata'};whichElectrodes];
%whichElectrodes = {'GB23','GB35','GB36','GB37', 'GB38' 'GB53', 'GB69','GB70','GB71','GB85','GB86','GB87'};
whichElectrodes = {'GB23', 'GB35', 'GB85', 'GB71'};

% Which stimulus conditions to plot? 
%whichTrials = {}; % if empty, plot average of all trials
whichTrials = {'CLENCH'};
%whichTrials = {'D', 'Y', 'V', 'F'};
%whichTrials = {'thumb','index','middle','ring','pinky'};

% Plot both broadband and evoked response per electrode
% 'out' contains the plotted time courses and SEs

% PLOT time course for each condition
% specs
specs = [];
specs.dataTypes          = {'broadband'};
specs.smoothingLevelInMs = [];
specs.collapseTrialTypes = 'no';
specs.baselineType       = 'selectedtrials';%'selectedtrials';
% plot specs
specs.plot.colorMap      = 'parula';
specs.plot.nSubPlots     = [];
specs.plot.addEccToTitle = 'yes';
specs.plot.showMax       = 'no';
specs.plot.fontSize      = 18;
specs.plot.Xlim          = [trials.time(1) trials.time(end)];


[out] = ecog_plotTimecourses(trials, whichElectrodes, whichTrials, specs);

%% TOP HALF oF HD GRID
v = 1:16:113;
inx = [v v+1 v+2 v+3 v+4 v+5 v+6 v+7];
[out] = ecog_plotTimecourses(trials, whichElectrodes(inx), whichTrials, specs);

%% BOTTOM HALF OF HD GRID
% inx = inx + 8; %inx = inx(1:length(inx)-2);
% [out] = ecog_plotTimecourses(trials, whichElectrodes(inx), whichTrials, specs);

%% OLD
%inx = [1:8 17:24 33:40 49:56 65:72 81:88 97:104 113:120];
%elecInx = 1:63:length(whichElectrodes)+63;
% for ii = 1:length(elecInx)-1
%     [out] = ecog_plotTimecourses(trials, whichElectrodes(elecInx(ii):elecInx(ii+1)-1), whichTrials, specs);
% end



