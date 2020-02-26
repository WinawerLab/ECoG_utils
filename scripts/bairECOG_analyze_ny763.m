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
%task              = 'sixcatlocisidiff';
%task              = 'sixcatloctemporal';
task              = 'prf';
specs.stim_names  = {'BLANK', 'VERTICAL', 'HORIZONTAL', 'DIAGONAL'};

%specs.stim_names  = {'TWOPULSE-1', 'TWOPULSE-2', 'TWOPULSE-3', 'TWOPULSE-4', 'TWOPULSE-5', 'TWOPULSE-6'};
%specs.stim_names  = {'TWOPULSE-1-0', 'TWOPULSE-2-0', 'TWOPULSE-3-0', 'TWOPULSE-4-0', 'TWOPULSE-5-0', 'TWOPULSE-6-0'};
%specs.stim_names  = {'TWOPULSE-1-1', 'TWOPULSE-2-1', 'TWOPULSE-3-1', 'TWOPULSE-4-1', 'TWOPULSE-5-1', 'TWOPULSE-6-1'};
%specs.stim_names  = {'ONEPULSE-1', 'ONEPULSE-2', 'ONEPULSE-3', 'ONEPULSE-4', 'ONEPULSE-5', 'ONEPULSE-6'};

%specs.stim_names  = {'TWOPULSE-1-FACES', 'TWOPULSE-2-FACES', 'TWOPULSE-3-FACES', 'TWOPULSE-4-FACES', 'TWOPULSE-5-FACES', 'TWOPULSE-6-FACES'};
%specs.stim_names  = {'ONEPULSE-1-FACES', 'ONEPULSE-2-FACES', 'ONEPULSE-3-FACES', 'ONEPULSE-4-FACES', 'ONEPULSE-5-FACES', 'ONEPULSE-6-FACES'};
specs.chan_names  = {'G1', 'G2', 'PT2', 'PT4'};
%specs.plot_ylim   = [-1 6];

bidsEcogPlotTrials(projectDir, subject, session, task, [], [], [], specs, 1);
%bidsEcogPlotTrials(projectDir, subject, session, task);


%% fit PRFs


% Load and epoch the data
recomputeFlag = true;
subjects      = {subject};
tasks         = {'prf'};
epochTime     = [-0.2 0.6];
saveStr       = 'prfdata';
saveDir       = fullfile(projectDir, 'derivatives', 'ECoGPRF', sprintf('sub-%s', subject), session);
[data] = tde_getData(recomputeFlag, subjects, [], tasks, [], epochTime, [], saveStr, saveDir);


% Compute the PRF timecourses
doPlots     = true;
plotSaveDir = fullfile(projectDir, 'derivatives', 'ECoGFigures', sprintf('sub-%s', subject), 'prfs');
[data]      = tde_computePRFtimecourses(data, [], [], doPlots, plotSaveDir);

% Load the stimulus apertures
stimName = fullfile(tdeRootPath, 'prf_apertures', 'bar_apertures.mat');
load(stimName, 'bar_apertures');

%% 2: Model fitting

% Fit the PRF time courses with analyzePRF
tr             = 1;
opt.hrf        = 1;
opt.maxpolydeg = 0;
opt.xvalmode   = 0; 
opt.display    = 'off';

doPlots = true;
[results] = tde_fitPRFs(data, bar_apertures, opt, doPlots, saveDir, [], plotSaveDir);
