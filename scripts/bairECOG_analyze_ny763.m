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
%task              = 'prf';
%specs.stim_names  = {'BLANK', 'VERTICAL', 'HORIZONTAL', 'DIAGONAL'};

specs.stim_names  = {'TWOPULSE-1', 'TWOPULSE-2', 'TWOPULSE-3', 'TWOPULSE-4', 'TWOPULSE-5', 'TWOPULSE-6'};
%specs.stim_names  = {'ONEPULSE-1', 'ONEPULSE-2', 'ONEPULSE-3', 'ONEPULSE-4', 'ONEPULSE-5', 'ONEPULSE-6'};
%specs.stim_names  = {'SAMETWOPULSE-1',  'SAMETWOPULSE-2', 'SAMETWOPULSE-3',  'SAMETWOPULSE-4', 'SAMETWOPULSE-5',  'SAMETWOPULSE-6'};
%specs.stim_names  = {'DIFFTWOPULSE-1',  'DIFFTWOPULSE-2', 'DIFFTWOPULSE-3',  'DIFFTWOPULSE-4', 'DIFFTWOPULSE-5',  'DIFFTWOPULSE-6'};

%specs.stim_names  = {'FACES',  'BODIES', 'OBJECTS', 'BUILDINGS','SCENES','SCRAMBLED'};
%specs.stim_names  = {'FACES-ONEPULSE','BODIES-ONEPULSE', 'OBJECTS-ONEPULSE', 'BUILDINGS-ONEPULSE','SCENES-ONEPULSE','SCRAMBLED-ONEPULSE'};
%specs.stim_names  = {'FACES-TWOPULSE-5','BODIES-TWOPULSE-5', 'OBJECTS-TWOPULSE-5', 'BUILDINGS-TWOPULSE-5','SCENES-TWOPULSE-5','SCRAMBLED-TWOPULSE-5'};
%specs.stim_names  = {'FACES-SAMETWOPULSE-5',  'BODIES-SAMETWOPULSE-5', 'OBJECTS-SAMETWOPULSE-5', 'BUILDINGS-SAMETWOPULSE-5','SCENES-SAMETWOPULSE-5','SCRAMBLED-SAMETWOPULSE-5'};
%specs.stim_names  = {'FACES-DIFFTWOPULSE-5',  'BODIES-DIFFTWOPULSE-5', 'OBJECTS-DIFFTWOPULSE-5', 'BUILDINGS-DIFFTWOPULSE-5','SCENES-DIFFTWOPULSE-5','SCRAMBLED-DIFFTWOPULSE-5'};
%specs.stim_names  = {'FACES-SAMETWOPULSE-1',  'FACES-SAMETWOPULSE-2', 'FACES-SAMETWOPULSE-3',  'FACES-SAMETWOPULSE-4', 'FACES-SAMETWOPULSE-5',  'FACES-SAMETWOPULSE-6'};
%specs.stim_names  = {'FACES-DIFFTWOPULSE-1',  'FACES-DIFFTWOPULSE-2', 'FACES-DIFFTWOPULSE-3',  'FACES-DIFFTWOPULSE-4', 'FACES-DIFFTWOPULSE-5',  'FACES-DIFFTWOPULSE-6'};

%specs.chan_names  = {'G01', 'G02', 'PT02', 'PT04'};
specs.chan_names  = {'G01', 'PT02'}; 
specs.plot_ylim   = [-1 6];
specs.plot_cmap   = 'copper';
specs.plot_includelegend   = 0;

task              = 'sixcatloctemporal';
specs.stim_names  = {'TWOPULSE-1', 'TWOPULSE-2', 'TWOPULSE-3', 'TWOPULSE-4', 'TWOPULSE-5', 'TWOPULSE-6'};
[out1] = bidsEcogPlotTrials(projectDir, subject, session, task, [], [], [], specs, 0);

task              = 'sixcatlocisidiff';
%specs.stim_names  = {'SAMETWOPULSE-1',  'SAMETWOPULSE-2', 'SAMETWOPULSE-3',  'SAMETWOPULSE-4', 'SAMETWOPULSE-5',  'SAMETWOPULSE-6'};
specs.stim_names  = {'DIFFTWOPULSE-1',  'DIFFTWOPULSE-2', 'DIFFTWOPULSE-3',  'DIFFTWOPULSE-4', 'DIFFTWOPULSE-5',  'DIFFTWOPULSE-6'};
[out2] = bidsEcogPlotTrials(projectDir, subject, session, task, [], [], [], specs, 0);

%% Plot Repetition Suppression

t_ind = out1{1}.t >0;

Arep = squeeze(sum(out1{1}.ts(t_ind,:,:),1));
Anonrep = squeeze(sum(out2{1}.ts(t_ind,:,:),1));

RS = (Anonrep-Arep)./ Anonrep;
%RS = (Anonrep-Arep);

figure;hold on
plot([1:6], RS(:,1),'k');
plot([1:6], RS(:,2),'--k');
scatter([1:6], RS(:,1), 250, out1{1}.colors, 'filled');
scatter([1:6], RS(:,2), 250, out1{1}.colors, 'filled', 'd');
set(gca, 'Xlim', [0 7]);

t_ind = out1{1}.t >0 & out1{1}.t < 0.25;
peakstim1 = squeeze(max(out1{1}.ts(t_ind,:,:),[],1));

t_ind = out1{1}.t >0 & out1{1}.t > 0.25;
peakstim2 = squeeze(max(out1{1}.ts(t_ind,:,:),[],1));

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

%%%  TO DO make separate plot function for PRF single trials and timeseries

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
chan_ind = ecog_matchChannels({'G01', 'G02', 'PT02', 'PT03'}, channels.name);
%chan_ind = [];

coloropt = 1;
figureName = sprintf('%s_prftimecoursefits', subject);
ecog_plotPRFtsfits(data2fit, stimulus, results, channels, chan_ind)
set(gcf, 'Position', get(0,'screensize'));
set(findall(gcf,'-property','FontSize'),'FontSize',14)
%saveas(gcf, fullfile(plotSaveDir, figureName), 'png'); close;

figureName = sprintf('%s_prfs', subject);
ecog_plotPRFs(results, stimulus, channels,chan_ind,[],coloropt)  
set(gcf, 'Position', get(0,'screensize'));
set(findall(gcf,'-property','FontSize'),'FontSize',14)
%saveas(gcf, fullfile(plotSaveDir, figureName), 'png'); close;


