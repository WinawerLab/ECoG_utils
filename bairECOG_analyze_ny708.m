
% SCRIPT DESCRIPTION %

%% [0] Define paths and dataset specs %%

% Dataset specs
projectName = 'visual';
sub_label   = 'som708'; 
ses_label   = {'nyuecog01','nyuecog02'};

% Input paths specs
dataPth     = '/Volumes/server/Projects/BAIR/Data/BIDS/';

%% Load preprocessed data

for ii = 1:length(ses_label)
    % Output paths specs
    procDir = fullfile(dataPth, projectName, 'derivatives', 'preprocessed', sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label{ii}));
    dataName = fullfile(procDir, sprintf('sub-%s_ses-%s_epoched.mat', sub_label, ses_label{ii}));
    all{ii} = load(dataName);
end    

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
trials = all{1}.trials;

% FROM GIO/DORA: Electrodes on seizure are OC12, OC13, OC21, OC22.
ElecsToExclude  = trials.channels.name(contains(trials.channels.status, 'bad'));

% Which electrodes to plot? (Each electrode gets a subplot)
%whichElectrodes = unique([trials.viselec.benson14_varea.elec_labels trials.viselec.wang15_mplbl.elec_labels]); 
%whichElectrodes = setdiff(whichElectrodes,ElecsToExclude);

% SELECTING ONLY V1, V2, V3, V3a, V3b:
%whichElectrodes = visualelectrodes.benson14_varea.elec_labels(strncmp(visualelectrodes.benson14_varea.area_labels, 'V',1));
%whichElectrodes = visualelectrodes.benson14_varea.elec_labels;
%whichElectrodes = ecog_sortElectrodesonVisualArea(whichElectrodes,visualelectrodes.wang15_mplbl);
whichElectrodes = trials.channels.name(contains(trials.channels.type, 'ecog'));
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
%whichTrials = {'CRF-5','ONEPULSE-5', 'TWOPULSE-5'}; 
whichTrials = {'HRF'};
%whichTrials = {'SPARSITY'};

% PLOT time course for each condition
specs = [];
specs.dataTypes          = {'broadband'};
specs.smoothingLevelInMs = [];
specs.collapseTrialTypes = 'no';
specs.baselineType       = 'all';%'selectedtrials';
specs.plot.colorMap      = 'parula';
specs.plot.nSubPlots     = [];
specs.plot.addEccToTitle = 'no';
specs.plot.showMax       = 'no';

% Plot both broadband and evoked response per electrode
% 'out' contains the plotted time courses and SEs
[out] = ecog_plotTimecourses(trials, whichElectrodes, whichTrials, specs);
%[out] = ecog_plotTrials(trials, whichElectrodes, whichTrials);
%%
elecInx = 1:(length(whichElectrodes)/2):length(whichElectrodes)+(length(whichElectrodes)/3);
for ii = 1:length(elecInx)-1
    [out] = ecog_plotTimecourses(trials, whichElectrodes(elecInx(ii):elecInx(ii+1)-1), whichTrials, specs);
end

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


data = trials.evoked;

data_epoch = zeros(size(data,1), size(data,3), size(data,2));

for ii = 1:size(data,3)
    data_epoch(:,ii,:) = data(:,:,ii);
end

stims = ones(1,size(data_epoch,2));
reg_erp = 0;
fft_w = [];
fft_t = 1:size(data_epoch,3);
fft_ov = [];
srate = trials.fsample;
%hdr.Fs,0,hdr.Fs,hdr.Fs

%[~,f] = pwelch(squeeze(data_epoch(1,1,:)),0,0,trials.fsample,trials.fsample);

[f,data_epoch_spectra,data_epoch] = ecog_spectra(data_epoch,stims,fft_w,fft_t,fft_ov,srate,reg_erp);


%% plot

whichElectrodes = find(contains(trials.channels.name, trials.viselec.benson14_varea.elec_labels));


whichTrials0 = contains(trials.events.trial_name, {'HRF'});
whichTrials1 = contains(trials.events.trial_name, {'GRATING', 'PLAID', 'CIRCULAR'});
whichTrials2 = contains(trials.events.trial_name, {'HOUSE', 'FACE', 'LETTER'});

figure;
for ee = 1:length(whichElectrodes)
    subplot(2,7,ee);hold on
    toplot = squeeze(mean(data_epoch_spectra(whichElectrodes(ee),whichTrials0,:),2));
    plot(f,log10(toplot),'k', 'LineWidth',2);
    toplot = squeeze(mean(data_epoch_spectra(whichElectrodes(ee),whichTrials1,:),2));
    plot(f,log10(toplot),'r', 'LineWidth',2);
    toplot = squeeze(mean(data_epoch_spectra(whichElectrodes(ee),whichTrials2,:),2));
    plot(f,log10(toplot),'b', 'LineWidth',2);
    if ee == 1
        legend({'HRF', 'ALLGRATINGS', 'OBJECTS'});
    end
    xlabel('frequency');
    ylabel('power');
    title(trials.viselec.benson14_varea.elec_labels(ee));
	set(gca, 'FontSize', 12);
end

whichTrials0 = contains(trials.events.trial_name, {'GRATING'});
whichTrials1 = contains(trials.events.trial_name, {'PLAID'});
whichTrials2 = contains(trials.events.trial_name, {'CIRCULAR'});

figure;
for ee = 1:length(whichElectrodes)
    subplot(2,7,ee);hold on
    toplot = squeeze(mean(data_epoch_spectra(whichElectrodes(ee),whichTrials0,:),2));
    plot(f,log10(toplot),'k', 'LineWidth',2);
    toplot = squeeze(mean(data_epoch_spectra(whichElectrodes(ee),whichTrials1,:),2));
    plot(f,log10(toplot),'r', 'LineWidth',2);
    toplot = squeeze(mean(data_epoch_spectra(whichElectrodes(ee),whichTrials2,:),2));
    plot(f,log10(toplot),'b', 'LineWidth',2);
    if ee == 1
        legend({'GRATING', 'PLAID', 'CIRCULAR'});
    end
    xlabel('frequency');
    ylabel('power');
    title(trials.viselec.benson14_varea.elec_labels(ee));
	set(gca, 'FontSize', 12);
end





