projectDir = '/Volumes/server/Projects/BAIR/Data/BIDS/visual'; 
subject = 'som765';

%% common average reference
bidsEcogRereference(projectDir, subject);

%% broadband 
bands = [[70 80]; [80 90]; [90 100]; [100 110]; [130 140]; [140 150]; [150 160]; [160 170]; [190 200]];
bidsEcogBroadband(projectDir, subject, [], [], [], bands);

%% plot
clear specs;

session           = 'nyuecog01';
task              = 'sixcatlocdiffisi';
%task              = 'sixcatloctemporal';
%task              = 'prf';
%specs.stim_names  = {'BLANK', 'VERTICAL', 'HORIZONTAL', 'DIAGONAL'};
%specs.stim_names  = {'FACES', 'BODIES', 'OBJECTS', 'BUILDINGS', 'SCENES', 'SCRAMBLED'};
%specs.stim_names  = {'ONEPULSE-1', 'ONEPULSE-2', 'ONEPULSE-3', 'ONEPULSE-4', 'ONEPULSE-5', 'ONEPULSE-6'};

specs.stim_names  = {'TWOPULSE-1', 'TWOPULSE-2', 'TWOPULSE-3', 'TWOPULSE-4', 'TWOPULSE-5', 'TWOPULSE-6'};
%specs.stim_names  = {'TWOPULSE-1-0', 'TWOPULSE-2-0', 'TWOPULSE-3-0', 'TWOPULSE-4-0', 'TWOPULSE-5-0', 'TWOPULSE-6-0'};
%specs.stim_names  = {'TWOPULSE-1-1', 'TWOPULSE-2-1', 'TWOPULSE-3-1', 'TWOPULSE-4-1', 'TWOPULSE-5-1', 'TWOPULSE-6-1'};%specs.stim_names  = {'ONEPULSE-1', 'ONEPULSE-2', 'ONEPULSE-3', 'ONEPULSE-4', 'ONEPULSE-5', 'ONEPULSE-6'};
specs.chan_names = {'PT'};
specs.plot_ylim = [-1 6];
bidsEcogPlotTrials(projectDir, subject, session, task, [], [], [], specs, 1);
%bidsEcogPlotTrials(projectDir, subject, session, task);


%% fit PRFs

% Load and epoch the data
dataDir       = fullfile(projectDir, 'derivatives', 'ECoGBroadband');
plotSaveDir   = fullfile(projectDir, 'derivatives', 'ECoGFigures', sprintf('sub-%s', subject), 'prfs');

session       = 'nyuecog01';
tasks         = 'prf';
description   = 'broadband';
epochTime     = [-0.2 0.6];
pRRFtime_win  = [0.05 0.55];

[data, channels, events] = bidsEcogGetPreprocData(dataDir, subject, session, tasks, [], description);
[channels]  = bair_addVisualAtlasNamesToChannelTable(channels,[]);
[epochs, t] = ecog_makeEpochs(data, events.onset, epochTime, channels.sampling_frequency(1));  

% Convert to PRF timeseries
ts = ecog_createPRFtimeseries(epochs,t,pRRFtime_win,unique(events.stim_file_index), 1);

% Load the stimulus apertures
stimName = fullfile(tdeRootPath, 'prf_apertures', 'bar_apertures.mat');
load(stimName, 'bar_apertures');
bar_apertures = imresize(bar_apertures, [100 100], 'nearest');

%% 2: Model fitting

% Fit the PRF time courses with analyzePRF
tr                  = 1;
opt.hrf             = 1;
opt.maxpolydeg      = 0;
opt.xvalmode        = 0; 
opt.display         = 'off';

% analyzePRFdog
opt.gaussianmode    = 'OG';         % please don't forget this option otherwise analyzePRFdog uses the DoG model
opt.prfmodel        = 'fixexpt';
opt.typicalexpt     = 0.05;
opt.forcebounds     = 1;

% Reformat data and apertues into cells inputs for analyzePRF
data2fit = []; stimulus = [];
for jj = 1:size(ts,3)
    data2fit{jj} = ts(:,:,jj);
    stimulus{jj} = double(bar_apertures);
end

results = analyzePRFdog(stimulus, data2fit, tr, opt);

%%
%load('/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/ECoGPRF/sub-som763/nyuecog01/sub-som763_prffits.mat');
%load('/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/ECoGPRF/sub-som763/nyuecog01/sub-som763_prffits_20200226T134533.mat');
%load('/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/ECoGPRF/sub-som763/nyuecog01/sub-som763_prffits_20200226T141627.mat');
%load('/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/ECoGPRF/sub-som763/nyuecog01/sub-som763_prffits_20200226T144748.mat')
%load('/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/ECoGPRF/sub-som763/nyuecog01/sub-som763_prffits_20200226T145556.mat');
%load('/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/ECoGPRF/sub-som763/nyuecog01/sub-som763_prffits_20200226T150544.mat')

%chan_ind = ecog_matchChannels({'G1', 'G2', 'G9', 'G17', 'AT4', 'PT2', 'PT3', 'PT4'}, results.channels.name);
%chan_ind = ecog_matchChannels(, channels.name);
chan_ind = find(contains(channels.group, 'strip'));
%chan_ind = [];

coloropt = 1;
figureName = sprintf('%s_prftimecoursefits', subject);
ecog_plotPRFtsfits(data2fit, stimulus, results, channels, chan_ind)
set(gcf, 'Position', get(0,'screensize'));
set(findall(gcf,'-property','FontSize'),'FontSize',14)
saveas(gcf, fullfile(plotSaveDir, figureName), 'png'); close;

figureName = sprintf('%s_prfs', subject);
ecog_plotPRFs(results, stimulus, channels,chan_ind,[],coloropt)  
set(gcf, 'Position', get(0,'screensize'));
set(findall(gcf,'-property','FontSize'),'FontSize',14)
saveas(gcf, fullfile(plotSaveDir, figureName), 'png'); close;


