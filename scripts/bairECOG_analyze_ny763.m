projectDir = '/Volumes/server/Projects/BAIR/Data/BIDS/visual'; 
subject = 'som763';

%% common average reference
bidsEcogRereference(projectDir, subject);

%% broadband 
bands = [[70 80]; [80 90]; [90 100]; [100 110]; [130 140]; [140 150]; [150 160]; [160 170]; [190 200]];
bidsEcogBroadband(projectDir, subject, [], [], [], bands);

%% plot
clear specs;

session           = 'nyuecog01';
task              = 'sixcatlocisidiff';
%task              = 'sixcatloctemporal';
%task              = 'prf';
%specs.stim_names  = {'BLANK', 'VERTICAL', 'HORIZONTAL', 'DIAGONAL'};

%specs.stim_names  = {'TWOPULSE-1', 'TWOPULSE-2', 'TWOPULSE-3', 'TWOPULSE-4', 'TWOPULSE-5', 'TWOPULSE-6'};
%specs.stim_names  = {'TWOPULSE-1-0', 'TWOPULSE-2-0', 'TWOPULSE-3-0', 'TWOPULSE-4-0', 'TWOPULSE-5-0', 'TWOPULSE-6-0'};
%specs.stim_names  = {'TWOPULSE-1-1', 'TWOPULSE-2-1', 'TWOPULSE-3-1', 'TWOPULSE-4-1', 'TWOPULSE-5-1', 'TWOPULSE-6-1'};
%specs.stim_names  = {'ONEPULSE-1', 'ONEPULSE-2', 'ONEPULSE-3', 'ONEPULSE-4', 'ONEPULSE-5', 'ONEPULSE-6'};

specs.stim_names  = {'FACES',  'BODIES', 'OBJECTS', 'BUILDINGS','SCENES','SCRAMBLED'};
specs.stim_names  = {'FACES-ONEPULSE',  'BODIES-ONEPULSE', 'OBJECTS-ONEPULSE', 'BUILDINGS-ONEPULSE','SCENES-ONEPULSE','SCRAMBLED-ONEPULSE'};
%specs.stim_names  = {'FACES-TWOPULSE-5',  'BODIES-TWOPULSE-5', 'OBJECTS-TWOPULSE-5', 'BUILDINGS-TWOPULSE-5','SCENES-TWOPULSE-5','SCRAMBLED-TWOPULSE-5'};
specs.stim_names  = {'FACES-SAMETWOPULSE-5',  'BODIES-SAMETWOPULSE-5', 'OBJECTS-SAMETWOPULSE-5', 'BUILDINGS-SAMETWOPULSE-5','SCENES-SAMETWOPULSE-5','SCRAMBLED-SAMETWOPULSE-5'};
specs.stim_names  = {'FACES-DIFFTWOPULSE-5',  'BODIES-DIFFTWOPULSE-5', 'OBJECTS-DIFFTWOPULSE-5', 'BUILDINGS-DIFFTWOPULSE-5','SCENES-DIFFTWOPULSE-5','SCRAMBLED-DIFFTWOPULSE-5'};
specs.stim_names  = {'DIFFTWOPULSE-1',  'SAMETWOPULSE-1'};
specs.stim_names  = {'SAMETWOPULSE-1',  'SAMETWOPULSE-2', 'SAMETWOPULSE-3',  'SAMETWOPULSE-4', 'SAMETWOPULSE-5',  'SAMETWOPULSE-6'};
specs.stim_names  = {'TWOPULSE-1', 'TWOPULSE-2', 'TWOPULSE-3', 'TWOPULSE-4', 'TWOPULSE-5', 'TWOPULSE-6'};
specs.stim_names  = {'FACES-SAMETWOPULSE-1',  'FACES-SAMETWOPULSE-2', 'FACES-SAMETWOPULSE-3',  'FACES-SAMETWOPULSE-4', 'FACES-SAMETWOPULSE-5',  'FACES-SAMETWOPULSE-6'};
%specs.stim_names  = {'FACES-TWOPULSE-1', 'FACES-TWOPULSE-2', 'FACES-TWOPULSE-3', 'FACES-TWOPULSE-4', 'FACES-TWOPULSE-5', 'FACES-TWOPULSE-6'};
specs.stim_names  = {'FACES-DIFFTWOPULSE-1',  'FACES-DIFFTWOPULSE-2', 'FACES-DIFFTWOPULSE-3',  'FACES-DIFFTWOPULSE-4', 'FACES-DIFFTWOPULSE-5',  'FACES-DIFFTWOPULSE-6'};
%specs.stim_names  = {'TWOPULSE-1-FACES', 'TWOPULSE-2-FACES', 'TWOPULSE-3-FACES', 'TWOPULSE-4-FACES', 'TWOPULSE-5-FACES', 'TWOPULSE-6-FACES'};
%specs.stim_names  = {'ONEPULSE-1-FACES', 'ONEPULSE-2-FACES', 'ONEPULSE-3-FACES', 'ONEPULSE-4-FACES', 'ONEPULSE-5-FACES', 'ONEPULSE-6-FACES'};
%specs.chan_names  = {'G01', 'G02', 'PT02', 'PT04'};
specs.chan_names  = {'G01'}; 
%specs.chan_names  = {'G1', 'PT2'};
specs.plot_ylim   = [-1 15];

bidsEcogPlotTrials(projectDir, subject, session, task, [], [], [], specs, 0);
%bidsEcogPlotTrials(projectDir, subject, session, task);


%% fit PRFs

% Load and epoch the data
dataDir       = fullfile(projectDir, 'derivatives', 'ECoGBroadband');
session       = 'nyuecog01';
tasks         = 'prf';
description   = 'broadband';
epochTime     = [-0.2 0.6];

[data, channels, events] = bidsEcogGetPreprocData(dataDir, subject, session, tasks, [], description);
[epochs, t] = ecog_makeEpochs(data, events.onset, epochTime, channels.sampling_frequency(1));  

[epochs]    = ecog_averageEpochs(epochs, events, stim_names); 

saveStr       = 'prfdata';
saveDir       = fullfile(projectDir, 'derivatives', 'ECoGPRF', sprintf('sub-%s', subject), session);


% Compute the PRF timecourses
doPlots     = false;
plotSaveDir = fullfile(projectDir, 'derivatives', 'ECoGFigures', sprintf('sub-%s', subject), 'prfs');
[data]      = tde_computePRFtimecourses(data, [], 1, doPlots, plotSaveDir);

% Load the stimulus apertures
stimName = fullfile(tdeRootPath, 'prf_apertures', 'bar_apertures.mat');
load(stimName, 'bar_apertures');

%% 2: Model fitting

% Fit the PRF time courses with analyzePRF
tr                  = 1;
opt.hrf             = 1;
opt.maxpolydeg      = 0;
opt.xvalmode        = 1; 
opt.display         = 'off';

% analyzePRFdog
opt.gaussianmode    = 'OG';         % please don't forget this option otherwise analyzePRFdog uses the DoG model
opt.prfmodel        = 'fixexpt';
opt.typicalexpt     = 0.05;
opt.forcebounds     = 1;

doPlots = true;

[results] = tde_fitPRFs(data, bar_apertures, opt, doPlots, saveDir, [], plotSaveDir);

%%
%load('/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/ECoGPRF/sub-som763/nyuecog01/sub-som763_prffits.mat');
%load('/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/ECoGPRF/sub-som763/nyuecog01/sub-som763_prffits_20200226T134533.mat');
%load('/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/ECoGPRF/sub-som763/nyuecog01/sub-som763_prffits_20200226T141627.mat');
%load('/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/ECoGPRF/sub-som763/nyuecog01/sub-som763_prffits_20200226T144748.mat')
%load('/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/ECoGPRF/sub-som763/nyuecog01/sub-som763_prffits_20200226T145556.mat');
%load('/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/ECoGPRF/sub-som763/nyuecog01/sub-som763_prffits_20200226T150544.mat')

%chan_ind = ecog_matchChannels({'G1', 'G2', 'G9', 'G17', 'AT4', 'PT2', 'PT3', 'PT4'}, results.channels.name);
chan_ind = ecog_matchChannels({'G1', 'G2', 'PT2', 'PT3'}, results.channels.name);
%chan_ind = [];

coloropt = 1;
figureName = sprintf('%s_prftimecoursefits', subject);
ecog_plotPRFtsfits(data2fit, stimulus, results, results.channels, chan_ind)
%set(gcf, 'Position', get(0,'screensize'));
set(findall(gcf,'-property','FontSize'),'FontSize',14)
%saveas(gcf, fullfile(plotSaveDir, figureName), 'png'); close;

figureName = sprintf('%s_prfs', subject);
ecog_plotPRFs(results, stimulus, results.channels,chan_ind,[],coloropt)  
%set(gcf, 'Position', get(0,'screensize'));
set(findall(gcf,'-property','FontSize'),'FontSize',14)
%saveas(gcf, fullfile(plotSaveDir, figureName), 'png'); close;


