
tbUse ECoG_utils;

%% Define paths and dataset

% Dataset specs
projectName = 'visual';
sub_label   = 'som726'; 
ses_label   = 'nyuecog02';

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

%% Extract broadband response to each PRF stimulus

% Select a subset of electrodes to analyze
elecIndex  = contains(trials.channels.name, 'GB'); % HD grid
channels   = trials.channels(elecIndex,:);

% Select the prf events
trialIndex = contains(trials.events.task_name, 'prf');
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
timeIndex  = [0.05 0.55];

% Compute average broadband response in time window
PRFbb_mean = squeeze(mean(PRFbb(:,trials.time>timeIndex(1) & trials.time<timeIndex(2),:),2)); 

% Reshape to separate the two runs
PRFbb_mean = reshape(PRFbb_mean,[size(PRFbb,1) size(PRFbb,3)/2 2]);

%% Look at the data

% Define a channel to plot:
chanToPlot = {'GB117'}; % 'GB117' 'GB118' 'GB119', 'GB116'' 

% Plot
chanIndex = ecog_matchChannels(chanToPlot, channels.name);
figure;hold on
plot(PRFbb_mean(chanIndex,:,1), 'b', 'LineWidth', 2); % first run
plot(PRFbb_mean(chanIndex,:,2), 'r', 'LineWidth', 2); % second run

% Add title, axes labels, legends etc
set(gca, 'XTick', 1:1:size(PRFbb_mean,2), 'XTickLabel', events.trial_name([1:size(PRFbb_mean,2)]), 'XTickLabelRotation', 90, 'FontSize',8);
title(sprintf('%s w:%s b:%s [ecc = %0.1f]', channels.name{chanIndex}, channels.wangarea{chanIndex}, channels.bensonarea{chanIndex}, channels.bensoneccen{chanIndex}),'FontSize',18);
set(gca, 'XLim', [0 size(PRFbb_mean,2)+1])
xlabel('PRF stimulus',  'FontSize',14);
ylabel('broadband power','FontSize',14);
legend({'PRF run 1', 'PRF run 2'}, 'FontSize',14);
set(gcf, 'Position', [60 300 2000 1000]);


%% Proceed to fitting....

load('/Users/winawerlab/matlab/toolboxes/BAIRstimuli/stimuli/bar_apertures.mat');
bar_apertures = imresize(bar_apertures, [100 100], 'nearest');

% Inputs to analyzePRF
stimulus = {bar_apertures,bar_apertures};
data = {PRFbb_mean(:,:,1),PRFbb_mean(:,:,2)};
tr = 1;

opt.hrf = 1;
opt.maxpolydeg = 0;

% Run analyzePRF
results = analyzePRF(stimulus,data,tr,opt);

%% Plot
ecc = []; ind = [];
for cc = 1:120
    tmp = channels.bensoneccen{cc};
    if ~ischar(tmp)
        ecc(cc) = tmp;
    else
        ecc(cc) = nan;
    end
    
end

% % height from center of screen to top of somMacbook
% screen_height = 1280;
% viewing_distance = 50;
% radius_in_cm  = screen_height/2; % radius in cm
% radius_in_deg = rad2deg(atan(radius_in_cm./viewing_distance)); % radius in degrees
% 
% numpix  = 100;%;ones(size(subj))*101; % number of pixels in stimulus description (one side) as input to fitprf
% pix2deg = radius_in_deg*2./numpix;
% xC      = (numpix+1)/2; % center location in pixels
% yC      = xC;            
% 
% prf_ecc = pix2deg.*results.ecc;
% % % convert pixels to degrees
% % x =  pix2deg.*(x0-xC);
% % y = -pix2deg.*(y0-yC);
% % s =  pix2deg.*s0;

% The stimulus is 100 pixels (in both height and weight), and this corresponds to
% 16 degrees of visual angle.  To convert from pixels to degreees, we multiply
% by 10/100.
cfactor =16.6/100;
prf_ecc = results.ecc*cfactor;

ind = results.R2>25 & results.ecc < 20;
corr(prf_ecc(ind), ecc(ind)','rows','complete')
figure;scatter(ecc(ind),prf_ecc(ind));
ylabel('eccentricity from prf fit');
xlabel('eccentricity from benson template');

ang = [];
for cc = 1:120
    tmp = channels.bensonangle{cc};
    if ~ischar(tmp)
        ang(cc) = tmp;
    else
        ang(cc) = nan;
    end
end
ang = ang + 90;%ang = mod(90-(-ang),360);
ind = results.R2>30 &results.ecc < 50;
corr(results.ang(ind), ang(ind)','rows','complete')
figure;scatter(ang(ind),results.ang(ind));
ylabel('angle from prf fit');
xlabel('angle from benson template');
axis([0 360 0 360])

idx = results.R2 > 30;
figure, 
subplot(1,2,1);polarplot(deg2rad(results.ang), prf_ecc, 'x')
hold on, polarplot(deg2rad(results.ang(idx)), prf_ecc(idx), 'x'); title('analyzePRF');
rlim([0 25])
subplot(1,2,2); polarplot(deg2rad(ang), ecc, 'x')
hold on, polarplot(deg2rad(ang(idx)), ecc(idx), 'x'); title('benson14');
rlim([0 25])

figure, scatter(squeeze(results.params(1,2,:)), squeeze(results.params(1,1,:)),'x')
hold on, scatter(squeeze(results.params(1,2,idx)), squeeze(results.params(1,1,idx)),'x')
axis([0 100 0 100])
xlabel('x'), ylabel('y');

%% Visualize the location of each voxel's pRF
figure; hold on;
set(gcf,'Units','points','Position',[100 100 400 400]);
cmap = jet(size(results.ang,1));
for p=1:size(results.ang,1)
  if results.R2(p)>40
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
axis([-10 10 -10 10]);
straightline(0,'h','k-');       % line indicating horizontal meridian
straightline(0,'v','k-');       % line indicating vertical meridian
axis square;
set(gca,'XTick',-10:2:10,'YTick',-10:2:10);
xlabel('X-position (deg)');
ylabel('Y-position (deg)');

%%
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
  %polymatrix{p} = projectionmatrix([]);
end


%%
% Which voxel should we inspect?  Let's inspect the second voxel.
vx = chanIndex;

% For each run, collect the data and the model fit.  We project out polynomials
% from both the data and the model fit.  This deals with the problem of
% slow trends in the data.
datats = {};
modelts = {};
for p=1:length(data)
  datats{p} =  polymatrix{p}*data{p}(vx,:)';
  modelts{p} = polymatrix{p}*modelfun(results.params(1,:,vx),stimulusPP{p});
end

% Visualize the results
figure; hold on;
set(gcf,'Units','points','Position',[100 100 1000 100]);
plot(cat(1,datats{:}),'r-');
plot(cat(1,modelts{:}),'b-');
straightline(224*(1:2)+.5,'v','g-');
xlabel('PRF stimulus');
ylabel('Broadband response');
ax = axis;
%axis([.5 1200+.5 ax(3:4)]);
title('Data and model fit');

%%

% %% Compute an analogous measure from the spectra?
% 
% specs = [];
% specs.window    = 100;%200;
% specs.ov        = 50; %100
% specs.reg_erp   = 0;
% 
% % trials
% specs.t         = [0.05 0.55]; 
% [spectra_trials] = ecog_computeTrialSpectra(trials, specs);
% 
% % blanks
% specs.t         = [-1 -0.5]; 
% [spectra_blanks] = ecog_computeTrialSpectra(blank_trials, specs);
% 
% % concatenate trials and blanks
% spectra           = spectra_trials;
% spectra.events    = [spectra_trials.events;spectra_blanks.events];
% spectra.pwrspctrm = cat(2,spectra_trials.pwrspctrm, spectra_blanks.pwrspctrm);
% 
% % extract spectra
% PRFspectra = spectra.pwrspctrm(elecIndex,trialIndex,:);
% 
