
% SCRIPT DESCRIPTION %

%% [0] Define paths and dataset specs %%

% Dataset specs
projectName = 'visual';
sub_label   = 'som661'; 
ses_label   = 'nyuecog01';

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

%% Plot broadband time courses for visual electrodes

saveFigure = 0;

% Which stimulus conditions to plot? 

%whichTrials = {'SCENES','FACES','LETTERS'};
%whichTrials = {'GRATING', 'PLAID','CIRCULAR'};
%whichTrials = {'CRF','ONEPULSE', 'TWOPULSE'}; 
whichTrials = {'PRF'};

%whichTrials = {'CRF-1','CRF-2', 'CRF-3','CRF-4', 'CRF-5'};
%whichTrials = {'ONEPULSE-1','ONEPULSE-2', 'ONEPULSE-3','ONEPULSE-4', 'ONEPULSE-5', 'ONEPULSE-6'};
%whichTrials = {'TWOPULSE-1','TWOPULSE-2', 'TWOPULSE-3','TWOPULSE-4', 'TWOPULSE-5', 'TWOPULSE-6'};
%whichTrials = {'SPARSITY-1','SPARSITY-2', 'CRF-5', 'SPARSITY-3','SPARSITY-4'};

% PLOT time course for each condition
specs = [];
specs.dataTypes          = {'broadband'};
specs.smoothingLevelInMs = [];
specs.collapseTrialTypes = 'no';
specs.baselineType       = 'selectedtrials';
specs.plot.colorMap      = 'parula';
specs.plot.nSubPlots     = [];
specs.plot.addEccToTitle = 'yes';
specs.plot.showMax       = 'no';
specs.plot.XLim          = [-0.2 1];
specs.plot.YLim          = [];%[-0.5 20];

% % Plot a subset of electrodes
whichElectrodes = trials.channels.name(contains(trials.channels.name, 'Oc'));
[trials_out] = ecog_plotTimecourses(trials, whichElectrodes, whichTrials, specs);
trialName = [whichTrials{:}]; 
plotName = 'bbtimecourse';

if saveFigure == 1
    saveLoc     = fullfile(figPth, projectName, sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label{whichSession}), 'figures');
    figureName  = sprintf('sub-%s_ses-%s_elecs-%s_%s_%s', sub_label, ses_label{whichSession}, inxNames{ee},trialName,plotName);
    printnice(gcf, [1 150], saveLoc, figureName);
end
    
%% Plot spectra

saveFigure   = 0;

whichTask = {'soc'};

%whichTrials = {'BLANK','HRF'};
%whichTrials = {'BLANK', 'SCENES','FACES','LETTERS'};
%whichTrials = {'BLANK','GRATING', 'PLAID','CIRCULAR'};
%whichTrials = {'BLANK','CRF','ONEPULSE', 'TWOPULSE'}; 
%whichTrials = {'BLANK','PRF'};

%whichTrials = {'BLANK','CRF-1','CRF-2', 'CRF-3','CRF-4', 'CRF-5'};
%whichTrials = {'BLANK','ONEPULSE-1','ONEPULSE-2', 'ONEPULSE-3','ONEPULSE-4', 'ONEPULSE-5', 'ONEPULSE-6'};
%whichTrials = {'BLANK','TWOPULSE-1','TWOPULSE-2', 'TWOPULSE-3','TWOPULSE-4', 'TWOPULSE-5', 'TWOPULSE-6'};
whichTrials = {'BLANK','SPARSITY-1','SPARSITY-2', 'CRF-5', 'SPARSITY-3','SPARSITY-4'};

specs = [];
specs.plot.colorMap      = 'parula';
specs.plot.nSubPlots     = [];
specs.plot.addEccToTitle = 'yes';
specs.plot.XLim          = [];
specs.plot.YLim          = [10^-3 10^3];%[];

% Plot a subset of electrodes
whichElectrodes = trials.channels.name(contains(trials.channels.name, 'OC'));
[spectra_out] = ecog_plotSpectra(spectra, whichElectrodes, whichTrials, whichTask, specs);
trialName = [whichTrials{:}]; 
plotName = 'bbtimecourse';

if saveFigure == 1
    saveLoc     = fullfile(figPth, projectName, sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label{whichSession}), 'figures');
    figureName  = sprintf('sub-%s_ses-%s_elecs-%s_%s_%s', sub_label, ses_label{whichSession}, inxNames{ee},trialName,plotName);
    printnice(gcf, [1 150], saveLoc, figureName);
end


%% OLD CODE FROM BEFORE BIDS %%

SNRepoch = [0.1 0.4]; % e.g. 80 to 400 ms
signaltoplot = 'evoked';
tasklabel = 'soc';

% load NY648_Elec2Node; % elec names from SOM txt file are not in same
% format as in data.label which is annoying; manually labeling electrodes
% with a node match in retinotopic atlas regions for now...
eltocolor1 = {'MO_01', 'MO_02', 'MO_03', 'MO_04', 'G_03', 'G_04', 'G_09', 'G_10', 'G_11', 'G_12', 'G_17', 'G_18', ...
             'DLPA_07', 'DLPA_08', 'DPSL_05'};
eltocolor2 = {'IO_01', 'IO_02', 'IO_03', 'IO_04'};

[~,t1] = min(abs(trials.time-SNRepoch(1)));
[~,t2] = min(abs(trials.time-SNRepoch(2)));
%temp = squeeze(sum(abs(trials.(signaltoplot).(tasklabel)(:,t1:t2,:)),2));
temp = squeeze(mean(trials.(signaltoplot).(tasklabel)(:,t1:t2,:),2));
SNR1 = abs(mean(temp,2));
SNR2 = SNR1./mean(std(trials.(signaltoplot).(tasklabel)(:,1:52,:),0,3),2); % divide by prestim std
SNR3 = SNR1./std(temp,0,2); % divide by std across trials

eltoplot = good_channels;
tempdata = data;
tempdata.label = tempdata.label(eltoplot);

% xaxislabels
for ee = 1:length(eltoplot)
    elnames{ee} = data.label{eltoplot(ee)}(4:end-4);
end

SNR = SNR2;
figure('Name', ['SNR per electrode: ' tasklabel]);hold on

el1 = ecog_matchchannels(eltocolor1,tempdata);
bar(1:length(eltoplot),SNR(eltoplot,:));
bar(el1,SNR(eltoplot(el1)),'FaceColor', 'r');
el2 = ecog_matchchannels(eltocolor2,tempdata);
bar(el2,SNR(eltoplot(el2)),'FaceColor', 'g');

legend({'non-matched electrodes', 'matched electrodes (wang atlas)', 'ventral'});
set(gca,'XTick', 1:126, 'XTickLabel', elnames, 'XTickLabelRotation', 90);
axis tight; ylabel(signaltoplot);
title(['interval ' num2str(SNRepoch(1)*1000) ' to ' num2str(SNRepoch(2)*1000) ' ms: abs(mean)/std prestim']);

%% average and compute an SNR measure per stim category

tasklabel = 'soc'; 
signaltoplot = 'broadband';
eltoplot = 'OC_01';

usebootstrap = 'no'; % yes or no

el = ecog_matchchannels(eltoplot,data);
temp = squeeze(trials.(signaltoplot).(tasklabel)(el,:,:)); 

% compute CATEGORY averages (combining all FACES, LETTERS, SCENES, and ONE_PULSE stimuli)
catList = [1:12 13 17 21 25 31:36 37];
nCat = length(catList);
ga = nan(nCat,length(trials.time)); 
ulim = nan(nCat, length(trials.time));  
llim = nan(nCat, length(trials.time)); 

for cCat = 1:nCat-1
    categoryInx = find(ismember(stimuli.(tasklabel), catList(cCat):catList(cCat+1)-1));
    disp(length(categoryInx));
    switch usebootstrap
        case 'yes'
            for tt = 1:length(trials.time)
                [bootstat] = bootstrp(100,@mean,squeeze(temp(tt,categoryInx)));
                ga(cCat,tt) = median(bootstat);
                llim(cCat,tt) = quantile(sort(bootstat), 0.05);
                ulim(cCat,tt) = quantile(sort(bootstat), 0.95);
            end
        otherwise
             ga(cCat,:) = mean(temp(:,categoryInx),2);
             llim(cCat,:) = ga(cCat,:) - std(temp(:,categoryInx),0,2)';
             ulim(cCat,:) = ga(cCat,:) + std(temp(:,categoryInx),0,2)';
    end 
end

% plot
figure('Name', [tasklabel ' ' signaltoplot ' ' data.label{el}]); hold on
lim = max(max(abs(temp)));
for cCat = 1:nCat-1
    subplot(ceil(sqrt(nCat)),ceil(sqrt(nCat)),cCat); hold on;

    h = ciplot(llim(cCat,:),ulim(cCat,:),trials.time,'k',0.1);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    plot(trials.time, ga(cCat,:), 'k', 'LineWidth',2);
    ylim = [-lim-0.2*lim lim+0.2*lim];
    set(gca, 'YLim', ylim);
    line([0 0], ylim,'LineStyle', ':', 'Color', 'k');
    line([-0.1 0.5], [0 0],'LineStyle', ':', 'Color', 'k');
    xlabel('time (s)');
    axis tight
    title(category_names{catList(cCat)})
end
