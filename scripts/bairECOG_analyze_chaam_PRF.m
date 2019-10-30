
%% Paths specs

% Input path
dataPth     = '/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/ECoGCAR';

% Dataset specs
subject     = 'chaam'; 
session     = 'umcuecog';
task        = {'prf'};
description = 'reref';

%% Get visual area matches for this subject

fprintf('Computing matches with visual atlases...\n');
specs = [];
specs.pID           = subject; 
specs.plotmesh      = 'none';
BIDSformatted       = 1;
visualelectrodes    = electrode_to_nearest_node(specs, BIDSformatted);
 
%% Get the data

fprintf('Reading in the data...\n');
[data, channels, events] = bidsEcogGetPreprocData(dataPth, subject, session, task, [], description);
srate = channels.sampling_frequency(1);

%% Data preprocessing

% % Shift onsets
fprintf('[%s] This is a umcu patient: shifting the onsets\n',mfilename);
shiftInSeconds = 0.062; % 62 ms
shiftInSamples = round(shiftInSeconds*srate); 
events.onset = events.onset + shiftInSeconds;
events.event_sample = events.event_sample + shiftInSamples; 

% Add visual area names (W and B) ecc, angle, sigma to channels table
[channels] = bair_addVisualAtlasNamesToChannelTable(channels,visualelectrodes);

% Epoch the data
epochTime = [-0.3 0.85];
fprintf('Computing epochs...\n');
[epochs, epoch_t] = ecog_makeEpochs(data, events.event_sample, epochTime, srate);  
%[epochs] = ecog_normalizeEpochs(epochs, epoch_t, [-0.2 0], 'subtractwithintrial');

%% %% check ERPs / normalization

%chan_ind = find(contains(channels.name, 'Oc18'));
%trial_ind = contains(events.trial_name, {'GRATING', 'CIRCULAR', 'PLAID', 'SPARSITY', 'CRF-5'});
%figure, hold on, plot(epoch_t,epochs(:,trial_ind,chan_ind)), plot(epoch_t,mean(epochs(:,trial_ind,chan_ind),2), 'k', 'LineWidth', 2);
%figure, hold on, plot(epoch_t,epochs_norm(:,trial_ind,chan_ind)), plot(epoch_t,mean(epochs_norm(:,trial_ind,chan_ind),2), 'k', 'LineWidth', 2);

% plot ERPs
chan_ind = find(contains(channels.name, 'Oc'));
%chan_ind = find(contains(channels.name, {'Oc19', 'Oc24', 'Oc16', 'sT1'}));
trial_ind = contains(events.task_name, 'prf') & ~contains(events.trial_name, 'BLANK');
blank_ind = contains(events.task_name, 'prf') & contains(events.trial_name, 'BLANK');

nChan = length(chan_ind); nSubPlot = ceil(sqrt(nChan));

figure;
for ii = 1:nChan
    subplot(nSubPlot, nSubPlot,ii); hold on
    plot(epoch_t, mean(epochs(:,trial_ind,chan_ind(ii)),2),'r', 'LineWidth', 2);
	plot(epoch_t, mean(epochs(:,blank_ind,chan_ind(ii)),2),'k', 'LineWidth', 2);
    ylimits = get(gca, 'YLim');
    l1 = line([0 0],ylimits, 'Color', 'k', 'LineStyle', ':');
    l2 = line(epochTime,[0 0],'Color', 'k', 'LineStyle', ':');
    title(channels.name(chan_ind(ii))); axis tight
end

%% Compute spectra using Welch

fft_w    = window(@hann,200);
fft_ov   = 100; % overlap
reg_erp  = 0;

% stim on
t    = [0 0.5]; 
fft_t = epoch_t > t(1) & epoch_t < t(2); 

% Dora's ecog_spectra function expects channels X epochs X time
data_epoch = permute(epochs, [3 2 1]);
stims = ones(1,size(data_epoch,2));

% compute spectra
[f,spectra] = ecog_spectra(data_epoch,stims,fft_w,fft_t,fft_ov,srate,reg_erp);

%% plot spectra
chan_ind = find(contains(channels.name, 'Oc'));
%chan_ind = find(contains(channels.name, {'Oc19', 'Oc24', 'Oc16', 'sT1'}));
trial_ind = contains(events.task_name, 'prf') & ~contains(events.trial_name, 'BLANK');
blank_ind = contains(events.task_name, 'prf') & contains(events.trial_name, 'BLANK');

nChan = length(chan_ind); nSubPlot = ceil(sqrt(nChan));

f_ind = 2:100;

figure;
for ii = 1:nChan
    subplot(nSubPlot, nSubPlot,ii); hold on
    chan_spectra = squeeze(spectra(chan_ind(ii),:,:));
    p1 = plot(f(f_ind), chan_spectra(trial_ind,f_ind), 'Color', [1 0.8 0.8], 'LineWidth', 1, 'LineStyle', ':');
    p2 = plot(f(f_ind), chan_spectra(blank_ind,f_ind), 'Color', [0.8 0.8 0.8], 'LineWidth', 1, 'LineStyle', ':');
    p1_m = plot(f(f_ind), mean(chan_spectra(trial_ind,f_ind)), 'r', 'LineWidth', 3);
    p2_m = plot(f(f_ind), mean(chan_spectra(blank_ind,f_ind)), 'k:', 'LineWidth', 3);
    %p1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    %p2.Annotation.LegendInformation.IconDisplayStyle = 'off';
    %legend({'PRF', 'BLANK'});
    set(gca, 'YScale', 'log', 'YLim', [10^-2.5 10^2.5]);
    title(channels.name(chan_ind(ii)));
end

figure;
for ii = 1:nChan
    subplot(nSubPlot, nSubPlot,ii); hold on
    chan_spectra = squeeze(spectra(chan_ind(ii),:,:));
    plot(f(f_ind), mean(chan_spectra(trial_ind,f_ind)) - mean(chan_spectra(blank_ind,f_ind)), 'k', 'LineWidth', 2);
    legend({'PRF-BLANK'});
    l1 = line([f(f_ind(1)) f(f_ind(end))],[0 0], 'Color', 'r', 'LineStyle', ':', 'LineWidth', 2);
    l1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    set(gca, 'YLim', [-10 10]);
    title(channels.name(chan_ind(ii)));
end

%% generate PRF time series

f_ind = [30:45 55:100];
prf_ts = geomean(spectra(:, :, f_ind),3);
nEvents = height(events);
prf_ts = reshape(prf_ts, [size(prf_ts,1) nEvents/2 2]);

chan_ind = find(contains(channels.name, 'Oc'));
%chan_ind = find(contains(channels.name, {'Oc19', 'Oc24', 'Oc16', 'sT1'}));

nChan = length(chan_ind); nSubPlot = ceil(sqrt(nChan));

figure;
for ii = 1:nChan
    subplot(nSubPlot, nSubPlot,ii); hold on
    plot(prf_ts(chan_ind(ii),:,1), 'r', 'LineWidth', 1)
    plot(prf_ts(chan_ind(ii),:,2), 'b', 'LineWidth', 2)
    if ii == 1
        legend({'run1', 'run2'});
    end
    axis tight
    title(channels.name(chan_ind(ii)));
end

%% fit PRF model

% Get bar apertures
load('/Users/winawerlab/matlab/toolboxes/BAIRstimuli/stimuli/bar_apertures.mat');
bar_apertures = imresize(bar_apertures, [100 100], 'nearest');

% Inputs to analyzePRF
stimulus = {bar_apertures};
data = {mean(prf_ts,3)}; % average over repeats
tr = 1;

opt.hrf = 1;
opt.maxpolydeg = 0;
opt.xvalmode = 0; % cross-validation across runs

% Run analyzePRF
results = analyzePRF(stimulus,data,tr,opt);

%% Plot R2
figure;plot(results.R2, 'LineWidth', 2);
set(gca, 'FontSize', 18, 'XTick', 1:8:height(channels));
xlabel('electrode inx'); ylabel('R2 (cross-validated)'); ylim([0 100]);

figure;histogram(results.R2, 100, 'FaceColor', 'b');
set(gca, 'XLim', [0 100], 'FontSize', 18);
xlabel('R2'); ylabel('number of elecs'); 

%% %% KK example code: Visualize the location of each voxel's pRF

% The stimulus is 100 pixels (in both height and weight), and this corresponds to
% 16.6 degrees of visual angle.  To convert from pixels to degreees, we multiply
% by 16.6/100.
cfactor = 16.6/100;

figure; hold on;
set(gcf,'Units','points','Position',[100 100 400 400]);
cmap = parula(nChan);
for p=1:nChan%size(results.ang,1)
  %if results.R2(p)>30
      xpos = results.ecc(chan_ind(p)) * cos(results.ang(chan_ind(p))/180*pi) * cfactor;
      ypos = results.ecc(chan_ind(p)) * sin(results.ang(chan_ind(p))/180*pi) * cfactor;
      ang = results.ang(chan_ind(p))/180*pi;
      sd = results.rfsize(chan_ind(p)) * cfactor;
      h = k_drawellipse(xpos,ypos,ang,2*sd,2*sd);  % circle at +/- 2 pRF sizes
      h.Annotation.LegendInformation.IconDisplayStyle = 'off';
      set(h,'Color',cmap(p,:),'LineWidth',2);
      set(scatter(xpos,ypos,'r.'),'CData',cmap(p,:));
  %end
end
k_drawrectangle(0,0,16.6,16.6,'k-');  % square indicating stimulus extent
axis([-20 20 -20 20]);
straightline(0,'h','k-');       % line indicating horizontal meridian
straightline(0,'v','k-');       % line indicating vertical meridian
axis square;
set(gca,'XTick',-20:2:20,'YTick',-20:2:20);
xlabel('X-position (deg)');
ylabel('Y-position (deg)');
legend(channels.name(chan_ind));
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

%% Comparison with benson

% The stimulus is 100 pixels (in both height and weight), and this corresponds to
% 16 degrees of visual angle.  To convert from pixels to degreees, we multiply
% by 10/100.
cfactor =16.6/100;
prf_ecc = results.ecc*cfactor;

ind = results.R2>30; %& results.ecc < 20;
corr(prf_ecc(ind), channels.bensoneccen(ind),'rows','complete')
figure;scatter(channels.bensoneccen(ind),prf_ecc(ind), 100, 'filled');
ylabel('eccentricity from analyzePRF');
xlabel('eccentricity from benson template');
set(gca, 'FontSize', 18)

%ang = ang + 90;%ang = mod(90-(-ang),360);
ind = results.R2>30; %&results.ecc < 50;
corr(results.ang(ind), channels.bensonangle(ind)+90,'rows','complete')
figure;scatter(channels.bensonangle(ind)+90,results.ang(ind),100,'filled');
ylabel('angle from analyzePRF');
xlabel('angle from benson template');
axis([0 360 0 360])
set(gca, 'FontSize', 18)

%% Plot results as polarplot, xy

idx = results.R2 > 50;
figure, 
subplot(1,2,1);polarplot(deg2rad(results.ang), prf_ecc, 'x')
hold on, polarplot(deg2rad(results.ang(idx)), prf_ecc(idx), 'x'); title('analyzePRF');
rlim([0 25]);
set(gca, 'FontSize', 18)

subplot(1,2,2); polarplot(deg2rad(channels.bensonangle+90), channels.bensoneccen, 'x')
hold on, polarplot(deg2rad(channels.bensonangle(idx)+90), channels.bensoneccen(idx), 'x'); title('benson14');
rlim([0 25])
set(gca, 'FontSize', 18)

figure, scatter(squeeze(results.params(1,2,:)), squeeze(results.params(1,1,:)),'x')
hold on, scatter(squeeze(results.params(1,2,idx)), squeeze(results.params(1,1,idx)),'x')
axis([0 100 0 100])
xlabel('x', 'FontSize', 28), ylabel('y', 'FontSize', 28);
set(gca, 'FontSize', 18)

%% Plot R2 in grid layout

results.R2 = resultsWBC.R2(:,1);
% v = 1:16:113;
% inx{1} = [v v+1 v+2 v+3 v+4 v+5 v+6 v+7]; % TOP HALF OF HD GRID
% inx{2} = inx{1} + 8; % BOTTOM HALF OF HD GRID
% test = [inx{1} inx{2}];

% R2 = reshape(results.R2, [16 8]); % this doesn't work because we have missing channels
fillIndex = setdiff(1:128,  [1 2 15 16 113 114 127 128]);

R2 = nan([16 8]); 
R2(fillIndex) = results.R2;
figure;imagesc(rot90(R2)); caxis([-10 80]);
c = colorbar; c.Label.String = 'R2';
set(gca, 'FontSize', 18, 'XTick', 1:1:16);
title('R2 grid layout - first half');

results.R2 = resultsWBC.R2(:,2);
% v = 1:16:113;
% inx{1} = [v v+1 v+2 v+3 v+4 v+5 v+6 v+7]; % TOP HALF OF HD GRID
% inx{2} = inx{1} + 8; % BOTTOM HALF OF HD GRID
% test = [inx{1} inx{2}];

% R2 = reshape(results.R2, [16 8]); % this doesn't work because we have missing channels
fillIndex = setdiff(1:128,  [1 2 15 16 113 114 127 128]);

R2 = nan([16 8]); 
R2(fillIndex) = results.R2;
figure;imagesc(rot90(R2)); caxis([-10 80]);
c = colorbar; c.Label.String = 'R2';
set(gca, 'FontSize', 18, 'XTick', 1:1:16);
title('R2 grid layout - second half');

