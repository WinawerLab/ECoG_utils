
% SCRIPT DESCRIPTION %

%% [0] Define paths and dataset specs %%

% Dataset specs
projectName = 'visual';
sub_label   = 'som718'; 
ses_label   = 'nyuecog01';

% Input paths specs
dataPth     = '/Volumes/server/Projects/BAIR/Data/BIDS/';

%% Load preprocessed data

% Output paths specs
procDir = fullfile(dataPth, projectName, 'derivatives', 'preprocessed', sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label));
dataName = fullfile(procDir, sprintf('sub-%s_ses-%s_epoched.mat', sub_label, ses_label));
load(dataName);

%% Do electrode matching
dataDir = fullfile(dataPth, projectName, sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label{1}), 'ieeg');

specs = [];
specs.pID           = sub_label;% patient ID number
specs.atlasNames    = {};
specs.patientPool   = 'BAIR';
specs.plotmesh      = 'left';
specs.plotlabel     = 'yes';
visualelectrodes    = electrode_to_nearest_node(specs, dataDir);

%% Plot responses for visual electrodes

% Which electrodes to plot? (Each electrode gets a subplot)
elInx = trials.channels.name(contains(trials.channels.status, 'good'));
%elInx = trials.channels.name(contains(trials.channels.type, {'ecog', 'eeg'}));
%elInx = trials.channels.name(contains(trials.channels.name, 'DPMT'));

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
whichTrials = {'CRF','ONEPULSE', 'TWOPULSE'}; 
%whichTrials = {'HRF'};
%whichTrials = {'SPARSITY'};

% PLOT time course for each condition
specs = [];
specs.dataTypes          = {'broadband'};
specs.smoothingLevelInMs = [];
specs.collapseTrialTypes = 'no';
specs.baselineType       = 'selectedtrials';
specs.plot.colorMap      = 'parula';
specs.plot.nSubPlots     = [];
specs.plot.addEccToTitle = 'no';
specs.plot.showMax       = 'no';

% Plot both broadband and evoked response per electrode
% 'out' contains the plotted time courses and SEs
[out] = ecog_plotTimecourses(trials, elInx, whichTrials, specs);
%[out] = ecog_plotTrials(trials, whichElectrodes, whichTrials);

%% SPECTRA
% Input:
% data_epoch: channels X epochs X time
%
% stims: if repress_erp==1, the ERP calculated per condition can be regressed out.
%   The condition is indicated in stims, a vector with one value for each
%   epoch. If choosing to regress erp out, make sure data are baseline
%   corrected first, so the ERP can be calculated correctly!
%
% Output:
% data_epoch_spectra is channels X epochs X frequencies
%
% Example:
% [f,data_epoch_spectra,data_epoch] = ...
%     ecog_spectra(data_epoch,stims,fft_w,fft_t,fft_ov,srate,reg_erp)
%


% dora settings
fft_w = window(@hann,200); % window width
fft_ov = 100; % overlap
% do not regress ERP here, because regressing out average evoked response with only a few trials can hurt 
reg_erp = 0; 

% TRIALS
data = trials.evoked;
data_epoch = zeros(size(data,1), size(data,3), size(data,2));
for ii = 1:size(data,3)
    data_epoch(:,ii,:) = data(:,:,ii);
end

stims = ones(1,size(data_epoch,2));
fft_t = trials.time >0.05 & trials.time < 0.55; 
srate = trials.fsample;

[f,trials_spectra] = ecog_spectra(data_epoch,stims,fft_w,fft_t,fft_ov,srate,reg_erp);

% BLANKS
data = blank_trials.evoked;
data_epoch = zeros(size(data,1), size(data,3), size(data,2));
for ii = 1:size(data,3)
    data_epoch(:,ii,:) = data(:,:,ii);
end

stims = ones(1,size(data_epoch,2));
fft_t = blank_trials.time > -0.7 & blank_trials.time < -0.2; 
srate = trials.fsample;

[f,blank_trials_spectra] = ecog_spectra(data_epoch,stims,fft_w,fft_t,fft_ov,srate,reg_erp);


%% plot
whichElectrodes = trials.channels.name(contains(trials.channels.status, 'good'));
elInx = ecog_matchChannels(whichElectrodes, trials);

whichTrials0 = contains(blank_trials.events.task_name, {'hrfpattern'});
whichTrials1 = contains(trials.events.trial_name, {'HRF'});
cmap = eval(specs.plot.colorMap);
colors = cmap(1:round((length(cmap)/3)):length(cmap),:);

figure;
for ee = 1:length(elInx)
    subplot(out.subplotdim(1),out.subplotdim(2),ee);hold on
    toplot = squeeze(mean(blank_trials_spectra(elInx(ee),whichTrials0,:),2));
    plot(f,toplot,'k', 'LineWidth',2);
    toplot = squeeze(mean(trials_spectra(elInx(ee),whichTrials1,:),2));
    %plot(f,toplot,'Color',colors(1,:), 'LineWidth',2);
    plot(f,toplot,'c', 'LineWidth',2);
    if ee == 1
        legend({'BLANK','HRF'});
        xlabel('frequency');
        ylabel('power');
    end
    title(whichElectrodes(ee));
	set(gca, 'FontSize', 10);
    set(gca, 'Xlim', [10 180], 'Yscale', 'log', 'YLim', [10^-2.5 10^2.5]);
end
set(gcf, 'Position', [150 100 2000 1250]);

%% plot
whichElectrodes = trials.channels.name(contains(trials.channels.status, 'good'));
elInx = ecog_matchChannels(whichElectrodes, trials);

whichTrials0 = contains(blank_trials.events.task_name, {'spatialobject'});
whichTrials1 = contains(trials.events.trial_name, {'HOUSES'});
whichTrials2 = contains(trials.events.trial_name, {'FACES'});
whichTrials3 = contains(trials.events.trial_name, {'LETTERS'});
cmap = eval(specs.plot.colorMap);
colors = cmap(1:round((length(cmap)/3)):length(cmap),:);

figure;
for ee = 1:length(elInx)
    subplot(out.subplotdim(1),out.subplotdim(2),ee);hold on
    toplot = squeeze(mean(blank_trials_spectra(elInx(ee),whichTrials0,:),2));
    plot(f,toplot,'k', 'LineWidth',2);
    toplot = squeeze(mean(trials_spectra(elInx(ee),whichTrials1,:),2));
    plot(f,toplot,'Color',colors(1,:), 'LineWidth',2);
    toplot = squeeze(mean(trials_spectra(elInx(ee),whichTrials2,:),2));
    plot(f,toplot,'Color',colors(2,:), 'LineWidth',2);
    toplot = squeeze(mean(trials_spectra(elInx(ee),whichTrials3,:),2));
    plot(f,toplot,'Color',colors(3,:), 'LineWidth',2);
    if ee == 1
        %legend({'PRF-BLANK','GRATING', 'PLAID', 'CIRCULAR'});
        %legend({'BLANK','HOUSESFACESLETTERS'});
        legend({'BLANK','HOUSES','FACES','LETTERS'});
        xlabel('frequency');
        ylabel('power');
    end
    title(whichElectrodes(ee));
	set(gca, 'FontSize', 10);
    set(gca, 'Xlim', [10 180], 'Yscale', 'log', 'YLim', [10^-2.5 10^2.5]);
end
set(gcf, 'Position', [150 100 2000 1250]);
%%
whichElectrodes = trials.channels.name(contains(trials.channels.status, 'good'));
elInx = ecog_matchChannels(whichElectrodes, trials);

whichTrials0 = contains(blank_trials.events.task_name, {'spatialpattern'});
whichTrials1 = contains(trials.events.trial_name, {'GRATING'});
whichTrials2 = contains(trials.events.trial_name, {'PLAID'});
whichTrials3 = contains(trials.events.trial_name, {'CIRCULAR'});
cmap = eval(specs.plot.colorMap);
colors = cmap(1:round((length(cmap)/3)):length(cmap),:);

figure;
for ee = 1:length(elInx)
    subplot(out.subplotdim(1),out.subplotdim(2),ee);hold on
    toplot = squeeze(mean(blank_trials_spectra(elInx(ee),whichTrials0,:),2));
    plot(f,toplot,'k', 'LineWidth',2);
    toplot = squeeze(mean(trials_spectra(elInx(ee),whichTrials1,:),2));
    plot(f,toplot,'Color',colors(1,:), 'LineWidth',2);
    toplot = squeeze(mean(trials_spectra(elInx(ee),whichTrials2,:),2));
    plot(f,toplot,'Color',colors(2,:), 'LineWidth',2);
    toplot = squeeze(mean(trials_spectra(elInx(ee),whichTrials3,:),2));
    plot(f,toplot,'Color',colors(3,:), 'LineWidth',2);
    if ee == 1
        %legend({'PRF-BLANK','GRATING', 'PLAID', 'CIRCULAR'});
        %legend({'BLANK','HOUSESFACESLETTERS'});
        legend({'BLANK','GRATING','PLAID','CIRCULAR'});
        xlabel('frequency');
        ylabel('power');
    end
    title(whichElectrodes(ee));
	set(gca, 'FontSize', 10);
    set(gca, 'Xlim', [10 180], 'Yscale', 'log', 'YLim', [10^-2.5 10^2.5]);
end
set(gcf, 'Position', [150 100 2000 1250]);
