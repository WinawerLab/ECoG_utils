
% SCRIPT DESCRIPTION %

%% [0] Define paths and dataset specs %%

% Dataset specs
projectName = 'visual';
sub_label   = 'som726'; 
ses_label   = {'nyuecog02','nyuecog03'};

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

whichSession = 1;
trials = all{whichSession}.trials;
blanks = all{whichSession}.blank_trials;

%% Compute spectra

specs = [];
specs.window    = 200;
specs.ov        = 100;
specs.reg_erp   = 0;

% trials
specs.t         = [0.05 0.55]; 
[spectra_trials] = ecog_computeTrialSpectra(trials, specs);

% blanks
specs.t         = [-1 -0.5]; 
[spectra_blanks] = ecog_computeTrialSpectra(blanks, specs);

% concatenate trials and blanks
spectra         = spectra_trials;
spectra.events  = [spectra_trials.events;spectra_blanks.events];
spectra.pwrspct = cat(2,spectra_trials.pwrspct, spectra_blanks.pwrspct);

%% Plot broadband time courses for visual electrodes

saveFigure = 0;

% Which stimulus conditions to plot? 

%whichTrials = {''}; % if empty, plot average of all trials

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
specs.plot.YLim          = []%[-0.5 20];

% % Plot a subset of electrodes
whichElectrodes = {'GB102', 'GB118', 'GB103', 'GB119'};
[out] = ecog_plotTimecourses(trials, whichElectrodes, whichTrials, specs);

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

plotName    = 'bbtimecourse_scaled';
trialName   = [whichTrials{:}]; % 'CRF' 'ONEPULSE';'TWOPULSE'; 'SPARSITY'

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

whichTask = {'hrfpattern'};

whichTrials = {'BLANK','HRF'};
%whichTrials = {'BLANK', 'HOUSES','FACES','LETTERS'};
%whichTrials = {'BLANK','GRATING', 'PLAID','CIRCULAR'};
%whichTrials = {'BLANK','CRF','ONEPULSE', 'TWOPULSE'}; 
%whichTrials = {'DIAGONAL'};

%whichTrials = {'CRF-1','CRF-2', 'CRF-3','CRF-4', 'CRF-5'};
%whichTrials = {'ONEPULSE-1','ONEPULSE-2', 'ONEPULSE-3','ONEPULSE-4', 'ONEPULSE-5', 'ONEPULSE-6'};
%whichTrials = {'TWOPULSE-1','TWOPULSE-2', 'TWOPULSE-3','TWOPULSE-4', 'TWOPULSE-5', 'TWOPULSE-6'};
%whichTrials = {'SPARSITY-1','SPARSITY-2', 'CRF-5', 'SPARSITY-3','SPARSITY-4'};

specs = [];
specs.plot.colorMap      = 'parula';
specs.plot.nSubPlots     = [];
specs.plot.addEccToTitle = 'yes';
specs.plot.XLim          = [];
specs.plot.YLim          = [10^-3 10^3];%[];

% Plot a subset of electrodes
whichElectrodes = {'GB102', 'GB118', 'GB103', 'GB119'};
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

for ee = 1:length(inx)

   [spectra_out{ee}] = ecog_plotSpectra(spectra, gridList(inx{ee}), whichTrials, whichTask, specs);

   if saveFigure == 1
        saveLoc     = fullfile(figPth, projectName, sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label{whichSession}), 'figures');
        figureName  = sprintf('sub-%s_ses-%s_elecs-%s_%s_%s', sub_label, ses_label{whichSession}, inxNames{ee},trialName,plotName);
        printnice(gcf, [1 150], saveLoc, figureName);
   end
end

