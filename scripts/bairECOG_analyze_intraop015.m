
% SCRIPT DESCRIPTION %

%% [0] Define paths and dataset specs %%

% Dataset specs
projectName = 'visual';
sub_label   = 'intraop015'; 
ses_label   = 'UMCUOR';

% Input paths 
projectDir  = '/Volumes/server/Projects/BAIR/';
dataPth     = fullfile(projectDir, 'Data', 'BIDS');

% Output path 
figPth      = fullfile(projectDir, 'Analyses');

%% Load preprocessed data

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
specs.t         = [-1.5 -1]; 
[spectra_blanks] = ecog_computeTrialSpectra(blank_trials, specs);

% concatenate trials and blanks
spectra           = spectra_trials;
spectra.events    = [spectra_trials.events;spectra_blanks.events];
spectra.pwrspctrm = cat(2,spectra_trials.pwrspctrm, spectra_blanks.pwrspctrm);

%% Plot broadband time courses for visual electrodes

saveFigure = 1;

% Which stimulus conditions to plot? 
whichTrials = {'ONEPULSE-1','ONEPULSE-2', 'ONEPULSE-3','ONEPULSE-4', 'ONEPULSE-5', 'ONEPULSE-6'};
%whichTrials = {'TWOPULSE-1','TWOPULSE-2', 'TWOPULSE-3','TWOPULSE-4', 'TWOPULSE-5', 'TWOPULSE-6'};
%whichTrials = {'VERTICAL','HORIZONTAL', 'DIAGONAL'};

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
specs.plot.XLim          = [-0.2 1.5];
specs.plot.YLim          = [-0.5 5];

% % Plot a subset of electrodes
whichElectrodes = {'chan97','chan98', 'chan99', 'chan116'};
ecog_plotTimecourses(trials, whichElectrodes, whichTrials, specs);

% % Plot the entire HD grid

% Which electrodes to plot? (Each electrode gets a subplot)
whichElectrodes = trials.channels.name(contains(trials.channels.name, 'chan'));
% Create a list of electrode names to match to the grid and which includes
% names for the electrodes of the grid that were not connected

% Plot response per electrode
v = 1:16:113;
inx{1} = [v v+1 v+2 v+3 v+4 v+5 v+6 v+7]; % TOP HALF OF HD GRID
inx{2} = inx{1} + 8; % BOTTOM HALF OF HD GRID
inxNames = {'HDtop', 'HDbottom'};

plotName    = 'bbtimecourse_notscaled';
trialName   = 'PRF';%%[whichTrials{:}]; % 'CRF' 'ONEPULSE';'TWOPULSE'; 'SPARSITY'

clear trials_out;
for ee = 1:length(inx)
    
    [trials_out{ee}] = ecog_plotTimecourses(trials, whichElectrodes(inx{ee}), whichTrials, specs);
    
    if saveFigure == 1
        saveLoc     = fullfile(figPth, projectName, sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label), 'figures');
        figureName  = sprintf('sub-%s_ses-%s_elecs-%s_%s_%s', sub_label, ses_label, inxNames{ee},trialName,plotName);
        printnice(gcf, [1 150], saveLoc, figureName);
    end
end

%% PRF

% Select a subset of electrodes to analyze
elecIndex  = contains(trials.channels.name, 'c'); % HD grid
channels   = trials.channels(elecIndex,:);

% Select the prf events
trialIndex = contains(trials.events.stim_file, 'prf');
events     = trials.events(trialIndex,:);

% Extract the prf data
PRFbb      = trials.broadband(elecIndex,:,trialIndex);

% Perform baseline correction

% % From ZHOU et al 2019: To convert the unit of the time-varying broadband
% to percent signal change in each electrode, we first averaged each
% broadband time series across epochs. We defined the first 200 ms prior to
% stimulus onset as the baseline period for the epoch-averaged time course,
% then we computed the percent signal change by dividing the entire 1200ms
% time course point-wise by the average of the baseline. To equalize the
% baseline across electrodes, we subtracted the baseline average from the
% entire time-course so each electrode has trial-averaged baseline 0.

base_range = (trials.time >= -0.2 & trials.time < 0);
m_base     = squeeze(median(mean(PRFbb(:, base_range, :), 2), 3));
PRFbb      = PRFbb./m_base-1;
    
% Define a time window over which to average the broadband timecourse
time_win  = [0.05 0.55];

% Compute average broadband response in time window
PRFbb_mean = squeeze(mean(PRFbb(:,trials.time>time_win(1) & trials.time<time_win(2),:),2)); 

% Reshape to separate the two runs
PRFbb_mean = reshape(PRFbb_mean,[size(PRFbb,1) size(PRFbb,3)/2 2]);

%% Look at the data

% Define a channel to plot:
chanToPlot = {'chan97','chan98', 'chan99', 'chan116'};

% Plot
chanIndex = ecog_matchChannels(chanToPlot, channels.name);
figure;hold on
colors = {'b','r'};
for ee = 1:length(chanToPlot{ii})
	chanIndex = ecog_matchChannels(chanToPlot{ee}, channels.name);
	subplot(2,2,ee); hold on;
    for ii = 1:size(PRFbb_mean,3)
        plot(PRFbb_mean(chanIndex,:,ii), colors{ii}, 'LineWidth', 2); % first run
    end

    % Add title, axes labels, legends etc
    %set(gca, 'XTick', 1:1:size(PRFbb_mean,2), 'XTickLabel', events.trial_name(1:size(PRFbb_mean,2)), 'XTickLabelRotation', 90, 'FontSize',8);
    %title(sprintf('%s w:%s b:%s [ecc = %0.1f]', channels.name{chanIndex}, channels.wangarea{chanIndex}, channels.bensonarea{chanIndex}, channels.bensoneccen(chanIndex)),'FontSize',28);
    title(channels.name{chanIndex},'FontSize',14);
    set(gca, 'XLim', [0 size(PRFbb_mean,2)+1], 'FontSize', 18)
    if ee == 1
        xlabel('PRF stimulus',  'FontSize',28);
        ylabel('broadband power','FontSize',28);
        legend({'PRF run 1', 'PRF run 2'}, 'FontSize',28);
    end
    set(gcf, 'Position', [150 100 2000 1250]);
end

%% across entire grid

% Plot response per electrode
v = 1:16:113;
inx{1} = [v v+1 v+2 v+3 v+4 v+5 v+6 v+7]; % TOP HALF OF HD GRID
inx{2} = inx{1} + 8; % BOTTOM HALF OF HD GRID
inxNames = {'HDtop', 'HDbottom'};

plotName    = 'prftimecourse';
trialName   = 'PRF';%%[whichTrials{:}]; % 'CRF' 'ONEPULSE';'TWOPULSE'; 'SPARSITY'
colors = {'b','r'};

for ii = 1:length(inx)
    
    figure;hold on
    
    for ee = 1:length(inx{ii})
        chanIndex = inx{ii}(ee);
        subplot(8,8,ee); hold on;
    
        for kk = 1:size(PRFbb_mean,3)
            plot(PRFbb_mean(chanIndex,:,kk), colors{kk}, 'LineWidth', 2); 
        end
    % Add title, axes labels, legends etc
    %set(gca, 'XTick', 1:1:size(PRFbb_mean,2), 'XTickLabel', events.trial_name(1:size(PRFbb_mean,2)), 'XTickLabelRotation', 90, 'FontSize',8);
        title(channels.name{chanIndex},'FontSize',14);
        set(gca, 'XLim', [0 size(PRFbb_mean,2)+1], 'FontSize', 14)
        set(gca, 'YLim', [0 100]);
        if ee == 1
            xlabel('PRF stimulus',  'FontSize',14);
            ylabel('sum of broadband power','FontSize',14);
            legend({'PRF run 1', 'PRF run 2'}, 'FontSize',14);
        end
    end
    set(gcf, 'Position', [150 100 2000 1250]);
    if saveFigure == 1
          saveLoc     = fullfile(figPth, projectName, sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label), 'figures');
          figureName  = sprintf('sub-%s_ses-%s_elecs-%s_%s_%s', sub_label, ses_label, inxNames{ii},trialName,plotName);
          printnice(gcf, [1 150], saveLoc, figureName);
    end
end

%% Proceed to fitting....

load('/Users/winawerlab/matlab/toolboxes/BAIRstimuli/stimuli/bar_apertures.mat');
bar_apertures = imresize(bar_apertures, [100 100], 'nearest');

% Inputs to analyzePRF
stimulus = {bar_apertures,bar_apertures};
data = {PRFbb_mean(:,:,1),PRFbb_mean(:,:,2)};
tr = 1;

opt.hrf = 1;
opt.maxpolydeg = 0;
opt.xvalmode = 0; 

% Run analyzePRF
results = analyzePRF(stimulus,data,tr,opt);


%% Plot R2
figure;plot(results.R2, 'LineWidth', 2);
set(gca, 'FontSize', 18, 'XTick', 1:8:height(channels));
xlabel('electrode inx'); ylabel('R2');

figure;histogram(results.R2, 100, 'FaceColor', 'b');
set(gca, 'XLim', [0 10], 'FontSize', 18);
xlabel('R2'); ylabel('number of elecs'); 

%% %% KK example code: Visualize the location of each voxel's pRF

% The stimulus is 100 pixels (in both height and weight), and this corresponds to
% 16.6 degrees of visual angle.  To convert from pixels to degreees, we multiply
% by 16.6/100.
cfactor = 16.6/100;

results.ang = resultsWOBC.ang(:,2);
results.ecc = resultsWOBC.ecc(:,2);
results.rfsize = resultsWOBC.rfsize(:,2);

figure; hold on;
set(gcf,'Units','points','Position',[100 100 400 400]);
cmap = jet(size(results.ang,1));
for p=1:size(results.ang,1)
  if results.R2(p)>30
      xpos = results.ecc(p) * cos(results.ang(p)/180*pi) * cfactor;
      ypos = results.ecc(p) * sin(results.ang(p)/180*pi) * cfactor;
      ang = results.ang(p)/180*pi;
      sd = results.rfsize(p) * cfactor;
      h = k_drawellipse(xpos,ypos,ang,2*sd,2*sd);  % circle at +/- 2 pRF sizes
      set(h,'Color',cmap(p,:),'LineWidth',2);
      set(scatter(xpos,ypos,'r.'),'CData',cmap(p,:));
  end
end
k_drawrectangle(0,0,16.6,16.6,'k-');  % square indicating stimulus extent
axis([-20 20 -20 20]);
straightline(0,'h','k-');       % line indicating horizontal meridian
straightline(0,'v','k-');       % line indicating vertical meridian
axis square;
set(gca,'XTick',-20:2:20,'YTick',-20:2:20);
xlabel('X-position (deg)');
ylabel('Y-position (deg)');

set(gca, 'FontSize', 18);
set(gcf, 'Position', [163   553   684   599]);

%% KK example code: plot time courses with fit

results = resultsWBC;
data = dataWBC;

% Define some variables
res = [100 100];                    % row x column resolution of the stimuli
resmx = 100;                        % maximum resolution (along any dimension)
hrf = results.options.hrf;          % HRF that was used in the model
degs = results.options.maxpolydeg;  % vector of maximum polynomial degrees used in the model

% Pre-compute cache for faster execution
[d,xx,yy] = makegaussian2d(resmx,2,2,2,2);

% Prepare the stimuli for use in the model
stimulusPP = {};
for p=1:length(stimulus)
  stimulusPP{p} = squish(stimulus{p},2)';  % this flattens the image so that the dimensionality is now frames x pixels
  stimulusPP{p} = [stimulusPP{p} p*ones(size(stimulusPP{p},1),1)];  % this adds a dummy column to indicate run breaks
end

% Define the model function.  This function takes parameters and stimuli as input and
% returns a predicted time-series as output.  Specifically, the variable <pp> is a vector
% of parameter values (1 x 5) and the variable <dd> is a matrix with the stimuli (frames x pixels).
% Although it looks complex, what the function does is pretty straightforward: construct a
% 2D Gaussian, crop it to <res>, compute the dot-product between the stimuli and the
% Gaussian, raise the result to an exponent, and then convolve the result with the HRF,
% taking care to not bleed over run boundaries.
modelfun = @(pp,dd) conv2run(posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res),makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*abs(pp(3))^2))); 0]) .^ posrect(pp(5)),hrf,dd(:,prod(res)+1));

% Construct projection matrices that fit and remove the polynomials.
% Note that a separate projection matrix is constructed for each run.
polymatrix = {};
for p=1:length(degs)
  polymatrix{p} = projectionmatrix(constructpolynomialmatrix(size(data{p},2),0:degs(p)));
end

% Pick a channel to inspect:
vx = chanIndex;

% % For each run, collect the data and the model fit.  We project out polynomials
% % from both the data and the model fit.  This deals with the problem of
% % slow trends in the data.
% datats = {};
% modelts = {};
% for p=1:length(data)
%   datats{p} =  polymatrix{p}*data{p}(vx,:)';
%   modelts{p} = polymatrix{p}*modelfun(results.params(1,:,vx),stimulusPP{p});
% end

% IRIS: do not project out polynomials:
for p=1:length(data)
    datats{p} = data{p}(vx,:)';
    modelts{p} = modelfun(results.params(1,:,vx),stimulusPP{p});
end

% Visualize the results
figure; hold on;
set(gcf,'Units','points','Position',[100 100 1000 100]);
plot(cat(1,datats{:}),'k-', 'LineWidth', 2);
plot(cat(1,modelts{:}),'r-','LineWidth', 2);
straightline(224*(1:2)+.5,'v','g-');
xlabel('PRF stimulus','FontSize', 28);
ylabel('Broadband response','FontSize', 28);
ax = axis;
%axis([.5 1200+.5 ax(3:4)]);
legend('data', 'model prediction');

set(gcf, 'Position', [60 300 2000 1000]);
set(gca, 'FontSize', 18)


%% Plot spectra

saveFigure   = 0;

whichTask = {'spatialobject'};

%whichTrials = {'BLANK','HRF'};
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

