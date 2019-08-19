
% SCRIPT DESCRIPTION %

%% [0] Define paths and dataset specs %%

% Dataset specs
projectName = 'visual';
sub_label   = 'som726'; 
ses_label   = {'nyuecog01','nyuecog02','nyuecog03'};

% Input paths 
projectDir  = '/Volumes/server/Projects/BAIR/';
dataPth     = fullfile(projectDir, 'Data', 'BIDS');

% Output path 
figPth      = fullfile(projectDir, 'Analyses');

%% Load preprocessed data

for ii = 1:length(ses_label)
    % Output paths specs
    procDir = fullfile(dataPth, projectName, 'derivatives', 'preprocessed', sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label{ii}));
    dataName = fullfile(procDir, sprintf('sub-%s_ses-%s_epoched.mat', sub_label, ses_label{ii}));
    all{ii}=load(dataName);
end    

%% Which session are we analyzing?

whichSession = 2;
trials = all{whichSession}.trials;
blanks = all{whichSession}.blank_trials;

%% Compute spectra

specs = [];
specs.window    = 100;%200;
specs.ov        = 50; %100
specs.reg_erp   = 0;

% trials
specs.t         = [0.05 0.55]; 
[spectra_trials] = ecog_computeTrialSpectra(trials, specs);

% blanks
specs.t         = [-1.5 -1]; 
[spectra_blanks] = ecog_computeTrialSpectra(blanks, specs);

% concatenate trials and blanks
spectra           = spectra_trials;
spectra.events    = [spectra_trials.events;spectra_blanks.events];
spectra.pwrspctrm = cat(2,spectra_trials.pwrspctrm, spectra_blanks.pwrspctrm);

%% Plot broadband time courses for visual electrodes

saveFigure = 0;

% Which stimulus conditions to plot? 

%whichTrials = {'HRF'};
whichTrials = {'HOUSES','FACES','LETTERS'};
%whichTrials = {'GRATING', 'PLAID','CIRCULAR'};
%whichTrials = {'CRF','ONEPULSE', 'TWOPULSE'}; 
%whichTrials = {'DIAGONAL'};

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
specs.plot.YLim          = [-0.5 20];

% % Plot a subset of electrodes
%whichElectrodes = {'GB102', 'GB118', 'GB103', 'GB119'};
%whichElectrodes = {'GB89'};
%ecog_plotTimecourses(trials, whichElectrodes, whichTrials, specs);

% % Plot the entire HD grid

% Which electrodes to plot? (Each electrode gets a subplot)
whichElectrodes = trials.channels.name(contains(trials.channels.name, 'GB'));
% Create a list of electrode names to match to the grid and which includes
% names for the electrodes of the grid that were not connected
gridList = [];
for ee = 1:128
    if ee < 10
        chanName = ['GB0' num2str(ee)];
    else
        chanName = ['GB' num2str(ee)];
    end
        
    if strmatch(chanName,whichElectrodes, 'exact')
        gridList{ee} = chanName;
    else
        gridList{ee} = [chanName '-nodata'];
    end
end

% Plot response per electrode
v = 1:16:113;
inx{1} = [v v+1 v+2 v+3 v+4 v+5 v+6 v+7]; % TOP HALF OF HD GRID
inx{2} = inx{1} + 8; % BOTTOM HALF OF HD GRID
inxNames = {'HDtop', 'HDbottom'};

plotName    = 'bbtimecourse_scaled';
trialName   = [whichTrials{:}]; % 'CRF' 'ONEPULSE';'TWOPULSE'; 'SPARSITY'

clear trials_out;
for ee = 1:length(inx)
    
    [trials_out{ee}] = ecog_plotTimecourses(trials, gridList(inx{ee}), whichTrials, specs);
    
    if saveFigure == 1
        saveLoc     = fullfile(figPth, projectName, sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label{whichSession}), 'figures');
        figureName  = sprintf('sub-%s_ses-%s_elecs-%s_%s_%s', sub_label, ses_label{whichSession}, inxNames{ee},trialName,plotName);
        printnice(gcf, [1 150], saveLoc, figureName);
    end
end

%% Plot spectra

saveFigure   = 0;

whichTask = {'spatialobject'};

%whichTrials = {'BLANK','HRF'};
whichTrials = {'BLANK', 'HOUSES','FACES','LETTERS'};
%whichTrials = {'BLANK','GRATING', 'PLAID','CIRCULAR'};
%whichTrials = {'BLANK','CRF','ONEPULSE', 'TWOPULSE'}; 
%whichTrials = {'BLANK','DIAGONAL', 'HORIZONTAL', 'VERTICAL'};

%whichTrials = {'BLANK','CRF-1','CRF-2', 'CRF-3','CRF-4', 'CRF-5'};
%whichTrials = {'BLANK','ONEPULSE-1','ONEPULSE-2', 'ONEPULSE-3','ONEPULSE-4', 'ONEPULSE-5', 'ONEPULSE-6'};
%whichTrials = {'BLANK','TWOPULSE-1','TWOPULSE-2', 'TWOPULSE-3','TWOPULSE-4', 'TWOPULSE-5', 'TWOPULSE-6'};
%whichTrials = {'SPARSITY-1','SPARSITY-2', 'CRF-5', 'SPARSITY-3','SPARSITY-4'};

specs = [];
specs.plot.colorMap      = 'parula';
specs.plot.nSubPlots     = [];
specs.plot.addEccToTitle = 'yes';
specs.plot.XLim          = [];
specs.plot.YLim          = [10^-3 10^3];%[];

% % Plot a subset of electrodes
whichElectrodes = {'GB89'};
[spectra_out] = ecog_plotSpectra(spectra, whichElectrodes, whichTrials, whichTask, specs);

%% Plot the entire HD grid

% Which electrodes to plot? (Each electrode gets a subplot)
whichElectrodes = trials.channels.name(contains(trials.channels.name, 'GB'));
% Create a list of electrode names to match to the grid and which includes
% names for the electrodes of the grid that were not connected
gridList = [];
for ee = 1:128
    if ee < 10
        chanName = ['GB0' num2str(ee)];
    else
        chanName = ['GB' num2str(ee)];
    end
        
    if strmatch(chanName,whichElectrodes, 'exact')
        gridList{ee} = chanName;
    else
        gridList{ee} = [chanName '-nodata'];
    end
end

% Plot response per electrode
v = 1:16:113;
inx{1} = [v v+1 v+2 v+3 v+4 v+5 v+6 v+7]; % TOP HALF OF HD GRID
inx{2} = inx{1} + 8; % BOTTOM HALF OF HD GRID
inxNames = {'HDtop', 'HDbottom'};

plotName     = 'spectra';
trialName    = [whichTrials{~contains(whichTrials, 'BLANK')}];

clear spectra_out;
for ee = 1:length(inx)

   [spectra_out{ee}] = ecog_plotSpectra(spectra, gridList(inx{ee}), whichTrials, whichTask, specs);

   if saveFigure == 1
        saveLoc     = fullfile(figPth, projectName, sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label{whichSession}), 'figures');
        figureName  = sprintf('sub-%s_ses-%s_elecs-%s_%s_%s', sub_label, ses_label{whichSession}, inxNames{ee},trialName,plotName);
        printnice(gcf, [1 150], saveLoc, figureName);
   end
end

%% compare eccentricity with face/house selectivity

whichElectrodes = trials.channels.name(contains(trials.channels.name, 'GB'));
timeInd         = [0.05 0.55]; 
trialNames      = {'FACES', 'HOUSES'};%, 'LETTERS'};

% Compute baseline
baseline_index = find(contains(trials.events.task_name, {'spatialobject'}));
baseline = mean(mean(trials.broadband(:,trials.time<0,baseline_index),2),3);

% Extract summed broadband response for each electrode and stimulus
BB_SUM = []; ecc = [];
for ee = 1:length(whichElectrodes)
    
    % Get electrode data
    elInx = ecog_matchChannels(whichElectrodes(ee),trials);
    elData = squeeze(trials.broadband(elInx,:,:));
    
    % Perform baseline correction
    elData = (elData - mean(elData(trials.time<0,:),1)) ./ baseline(elInx);
    
    % Extract response time course to each category
    BB = [];
    for ii = 1:length(trialNames)
        BB(:,ii) = median(elData(:, contains(trials.events.trial_name, trialNames{ii})),2);
    end
    
    % Compute sum over timeInd
    temp = sum(BB(trials.time>timeInd(1) & trials.time<timeInd(2),:));
    if any(temp>0)
        BB_SUM(ee,:) = temp;
    else
        BB_SUM(ee,:) = nan([1 length(trialNames)]);
    end
    
    ecc(ee) = trials.channels.bensoneccen(elInx);
end

% diffBB_SUM = [];
% diffBB_SUM(:,1) = BB_SUM(:,1) - mean(BB_SUM(:,[2 3]),2); % FACES - [HOUSES LETTERS]
% diffBB_SUM(:,2) = BB_SUM(:,2) - mean(BB_SUM(:,[1 3]),2); % HOUSES - [FACES LETTERS]
% diffBB_SUM(:,3) = BB_SUM(:,3) - mean(BB_SUM(:,[1 2]),2); % HOUSES - [FACES LETTERS]

% Compute selectivity
diffBB_SUM = [];
diffBB_SUM(:,1) = BB_SUM(:,1) - BB_SUM(:,2); % [FACES -  HOUSES]
diffBB_SUM(:,2) = BB_SUM(:,2) - BB_SUM(:,1); % [HOUSES - FACES]

%% compare eccentricity with face/house selectivity

% Bin the eccentricity values
edges = [0:2:10 32];
%edges = [0 1 2 3 4 6 8 10 32];
Y = discretize(ecc,edges);

xnames = [];
for ii = 1:length(unique(Y(~isnan(Y))))
    xnames{ii} = sprintf('%d - %d', edges(ii),edges(ii+1));
end

% Plot
figure;hold on

subplot(3,2,1); hold on;
scatter(ecc, BB_SUM(:,1),'b', 'filled')
scatter(ecc, BB_SUM(:,2),'r', 'filled')
%scatter(ecc, BB_SUM(:,3),'g', 'filled')
legend(trialNames);
ylabel('sum of mean broadband timecourse');
xlabel('eccentricity');
%set(gca, 'XTick', 1:max(Y), 'XTickLabel', xnames, 'XLim', [0 max(Y)+1]);

subplot(3,2,3); hold on;
scatter(Y, BB_SUM(:,1),'b', 'filled')
scatter(Y, BB_SUM(:,2),'r', 'filled')
%scatter(Y, BB_SUM(:,3),'g', 'filled')
legend(trialNames);
ylabel('summed broadband');
xlabel('eccentricity bin');
set(gca, 'XTick', 1:length(xnames), 'XTickLabel', xnames, 'XLim', [0  length(xnames)+1]);

subplot(3,2,2); hold on;
scatter(ecc, diffBB_SUM(:,1),'k', 'filled')
%scatter(ecc, diffBB_SUM(:,2),'r', 'filled')
%scatter(ecc, BB_SUM(:,3),'g', 'filled')
legend({'FACES-HOUSES'});
ylabel('difference in broadband');
xlabel('eccentricity');
%set(gca, 'XTick', 1:max(Y), 'XTickLabel', xnames, 'XLim', [0 max(Y)+1]);

subplot(3,2,4); hold on;
scatter(Y, diffBB_SUM(:,1),'k', 'filled')
%scatter(Y, diffBB_SUM(:,2),'r', 'filled')
%scatter(Y, BB_SUM(:,3),'g', 'filled')
legend({'FACES-HOUSES'});
ylabel('difference in broadband');
xlabel('eccentricity bin');
set(gca, 'XTick', 1:length(xnames), 'XTickLabel', xnames, 'XLim', [0  length(xnames)+1]);

% Compute and plot mean per bin
meanBB = []; seBB = [];
for ii = 1:max(Y)
    meanBB(ii,:) = nanmedian(diffBB_SUM(Y==ii,:));
    seBB(ii,:) = nanstd(diffBB_SUM(Y==ii,:), 0, 1)/sqrt(length(find(Y==ii)));
end

subplot(3,2,5); hold on;
p1 = plot(1:max(Y), meanBB(:,1),'k', 'LineWidth', 2);
%p2 = errorbar(1:max(Y), meanBB(:,1), seBB(:,1), 'k');
% for ii = 1:length(unique(Y))
%     nElecs = length(find(Y==ii));
%     p = plot(ones([nElecs 1])*ii,BB_SUM(Y==ii,1),'b.','MarkerSize', 30,'LineStyle', 'none');
%     p.Annotation.LegendInformation.IconDisplayStyle = 'off';
% end
legend('FACES-HOUSES');
ylabel('median of difference');
xlabel('eccentricity bin');
set(gca, 'XTick', 1:length(xnames), 'XTickLabel', xnames, 'XLim', [0  length(xnames)+1]);

% Boxplot
subplot(3,2,6); hold on;
boxplot(diffBB_SUM(:,1),Y);
ylabel('median of difference');
xlabel('eccentricity bin');
set(gca, 'XTick', 1:length(xnames), 'XTickLabel', xnames, 'XLim', [0  length(xnames)+1]);

set(gcf, 'Position', [400 400 1300 1000])
