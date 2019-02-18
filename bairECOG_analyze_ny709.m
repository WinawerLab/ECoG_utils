
% SCRIPT DESCRIPTION %

%% [0] Define paths and dataset specs %%

% Dataset specs
projectName = 'visual';
sub_label   = 'som709'; 
ses_label   = 'nyuecog01';

% Input paths specs
dataPth     = '/Volumes/server/Projects/BAIR/Data/BIDS/';

%% Load preprocessed data

% Output paths specs
procDir = fullfile(dataPth, projectName, 'derivatives', 'preprocessed', sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label));
dataName = fullfile(procDir, sprintf('sub-%s_ses-%s_epoched.mat', sub_label, ses_label));
load(dataName);

%% Do electrode matching
dataDir = fullfile(dataPth, projectName, sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label), 'ieeg');

specs = [];
specs.pID           = sub_label;% patient ID number
specs.atlasNames    = {};
specs.patientPool   = 'BAIR';
specs.plotmesh      = 'left';
specs.plotlabel     = 'yes';
visualelectrodes    = electrode_to_nearest_node(specs, dataDir);

%% Plot responses for visual electrodes

% FROM GIO/DORA: Electrodes on seizure are OC12, OC13, OC21, OC22.
ElecsToExclude  = trials.channels.name(contains(trials.channels.status, 'bad'));

% Which electrodes to plot? (Each electrode gets a subplot)
%whichElectrodes = unique([trials.viselec.benson14_varea.elec_labels trials.viselec.wang15_mplbl.elec_labels]); 
%whichElectrodes = setdiff(whichElectrodes,ElecsToExclude);

% SELECTING ONLY V1, V2, V3, V3a, V3b:
%whichElectrodes = visualelectrodes.benson14_varea.elec_labels(strncmp(visualelectrodes.benson14_varea.area_labels, 'V',1));
%whichElectrodes = visualelectrodes.benson14_varea.elec_labels;
%whichElectrodes = ecog_sortElectrodesonVisualArea(whichElectrodes,visualelectrodes.wang15_mplbl);
whichElectrodes = trials.channels.name(contains(trials.channels.type, {'seeg', 'ecog'}));
%whichElectrodes = {'O02', 'O03', 'O04', 'P03','P04','P05'};   %'O01','O05', 

% Which stimulus conditions to plot? 
%whichTrials = {''}; % if empty, plot average of all trials
%whichTrials = {'CRF','SPARSITY','ONEPULSE', 'TWOPULSE', 'GRATING', 'PLAID', 'CIRCULAR', 'HOUSES', 'FACES', 'LETTERS'}; % everything except prf
%whichTrials = {'CRF-1','CRF-2', 'CRF-3','CRF-4', 'CRF-5'};
%whichTrials = {'ONEPULSE-1','ONEPULSE-2', 'ONEPULSE-3','ONEPULSE-4', 'ONEPULSE-5', 'ONEPULSE-6'};
%whichTrials = {'TWOPULSE-1','TWOPULSE-2', 'TWOPULSE-3','TWOPULSE-4', 'TWOPULSE-5', 'TWOPULSE-6'};
%whichTrials = {'HOUSES','FACES', 'LETTERS'};
%whichTrials = {'GRATING', 'PLAID','CIRCULAR'};
%whichTrials = {'CRF', 'SPARSITY','ONEPULSE','GRATING', 'PLAID','CIRCULAR'};
%whichTrials = {'SPARSITY-1','SPARSITY-2', 'SPARSITY-3','SPARSITY-4'};
whichTrials = {'HRF','GRATING','PLAID', 'CIRCULAR','CRF-5','ONEPULSE-5', 'HOUSES', 'FACES', 'LETTERS'}; 
%whichTrials = {'HRF'};
%whichTrials = {'SPARSITY'};

% PLOT time course for each condition
specs = [];
specs.dataTypes          = {'broadband'};
specs.smoothingLevelInMs = [];
specs.collapseTrialTypes = 'yes';
specs.baselineType       = 'selectedtrials';%'selectedtrials';
specs.plot.colorMap      = 'parula';
specs.plot.nSubPlots     = [];
specs.plot.addEccToTitle = 'no';
specs.plot.showMax       = 'no';

% Plot both broadband and evoked response per electrode
% 'out' contains the plotted time courses and SEs
%[out] = ecog_plotTimecourses(trials, whichElectrodes, whichTrials, specs);
%[out] = ecog_plotTrials(trials, whichElectrodes, whichTrials);

%
elecInx = 1:(length(whichElectrodes)/2):length(whichElectrodes)+(length(whichElectrodes)/2);
for ii = 1:length(elecInx)-1
    [out] = ecog_plotTimecourses(trials, whichElectrodes(elecInx(ii):elecInx(ii+1)-1), whichTrials, specs);
end





