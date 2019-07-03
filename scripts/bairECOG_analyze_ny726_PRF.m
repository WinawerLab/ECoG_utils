
tbUse ECoG_utils;

%% [0] Define paths and dataset specs %%

% Dataset specs
projectName = 'visual';
sub_label   = 'som726'; 
ses_label   = 'nyuecog03';

% Input paths 
projectDir  = '/Volumes/server/Projects/BAIR/';
dataPth     = fullfile(projectDir, 'Data', 'BIDS');

% Output path 
figPth      = fullfile(projectDir, 'Analyses');

%% Load preprocessed data

% Output paths specs
procDir = fullfile(dataPth, projectName, 'derivatives', 'preprocessed', sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label));
dataName = fullfile(procDir, sprintf('sub-%s_ses-%s_epoched.mat', sub_label, ses_label));
load(dataName);

%% Compute spectra

specs = [];
specs.window    = 100;%200;
specs.ov        = 50; %100
specs.reg_erp   = 0;

% trials
specs.t         = [0.05 0.55]; 
[spectra_trials] = ecog_computeTrialSpectra(trials, specs);

% blanks
specs.t         = [-1 -0.5]; 
[spectra_blanks] = ecog_computeTrialSpectra(blank_trials, specs);

% concatenate trials and blanks
spectra           = spectra_trials;
spectra.events    = [spectra_trials.events;spectra_blanks.events];
spectra.pwrspctrm = cat(2,spectra_trials.pwrspctrm, spectra_blanks.pwrspctrm);

%% Extract broadband response to each PRF stimulus

% Select a subset of electrodes to analyze
elecIndex  = contains(trials.channels.name, 'GB'); % HD grid
channels   = trials.channels(elecIndex,:);

% Select the prf trials
trialIndex = contains(trials.events.task_name, 'prf');
events     = trials.events(trialIndex,:);

% Select the prf trials
PRFbb      = trials.broadband(elecIndex,:,trialIndex);
PRFspectra = spectra.pwrspctrm(elecIndex,trialIndex,:);

% Define a time window over which to average the broadband timecourse
timeIndex  = [0.05 0.55];

% Compute average broadband response in time window
PRFbb_mean = squeeze(mean(PRFbb(:,trials.time>timeIndex(1) & trials.time<timeIndex(2),:),2)); % should we take sum?

% Reshape to separate the two runs
PRFbb_mean = reshape(PRFbb_mean,[size(PRFbb,1) size(PRFbb,3)/2 2]);

% Define a channel to plot:
chanToPlot = {'GB117'}; % 'GB117' 'GB118' 'GB119', 'GB116'' 

% Plot 
chanIndex = ecog_matchChannels(chanToPlot, channels.name);
figure;hold on
plot(PRFbb_mean(chanIndex,:,1), 'b', 'LineWidth', 2); % first run
plot(PRFbb_mean(chanIndex,:,2), 'r', 'LineWidth', 2); % second run

% Add title, axes labels, legends etc
set(gca, 'XTick', 1:1:size(PRFbb_mean,2), 'XTickLabel', events.trial_name([1:size(PRFbb_mean,2)]), 'XTickLabelRotation', 90, 'FontSize',8);
title(sprintf('%s w:%s b:%s [ecc = %0.1f]', channels.name{chanIndex}, channels.wangarea{chanIndex}, channels.bensonarea{chanIndex}, channels.bensoneccen{chanIndex}),'FontSize',18);
set(gca, 'XLim', [0 size(PRFbb_mean,2)+1])
xlabel('PRF stimulus',  'FontSize',14);
ylabel('broadband power','FontSize',14);
legend({'PRF run 1', 'PRF run 2'}, 'FontSize',14);
set(gcf, 'Position', [60 300 2000 1000]);

% Compute an analogous measure from the spectra?

% Proceed to fitting....
% tbUse ECoG_prfanalysis




