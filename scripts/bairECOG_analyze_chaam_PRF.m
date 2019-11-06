
%% Paths specs

% Input path
dataPth     = '/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/ECoGCAR';
savePth     = '/Volumes/server/Projects/BAIR/Analyses/visual/sub-chaam/analyzePRF_nyu/';

% Dataset specs
subject     = 'chaam'; 
session     = 'umcuecog';
task        = {'prf'};
description = 'reref';

% Get bar apertures (necessary for PRF fitting)
%a = load('/Volumes/GoogleDrive/My Drive/temp/bar_apertures.mat');
a = load('/Users/winawerlab/matlab/toolboxes/BAIRstimuli/stimuli/bar_apertures.mat');

% Channels of interest
COI = {'Oc17', 'Oc18', 'Oc24', 'Oc16', 'sT1', 'Oc31'}; 

%% Get the data
fprintf('Reading in the data...\n');
[data, channels, events] = bidsEcogGetPreprocData(dataPth, subject, session, task, [], description);
srate = channels.sampling_frequency(1);
chan_inx = find(contains(channels.name, COI));

%% Get visual area matches for this subject

fprintf('Computing matches with visual atlases...\n');
specs = [];
specs.pID           = subject; 
specs.plotmesh      = 'none';
specs.plotlabel     = 'yes'; %specs.atlasNames    = {'benson14_varea'};
BIDSformatted       = 1;
visualelectrodes    = electrode_to_nearest_node(specs, BIDSformatted);

% Add visual area names (W and B) ecc, angle, sigma to channels table
[channels] = bair_addVisualAtlasNamesToChannelTable(channels,visualelectrodes);

%% Data preprocessing

%% Shift onsets (estimated delay with respect to NYU data)
fprintf('[%s] This is a umcu patient: shifting the onsets\n',mfilename);
shiftInSeconds = 0.062; % 62 ms
shiftInSamples = round(shiftInSeconds*srate); 
newevents = events;
newevents.onset = events.onset + shiftInSeconds;
newevents.event_sample = events.event_sample + shiftInSamples; 

%% Epoch the data
epochTime = [-0.3 0.85];
fprintf('Computing epochs...\n');
[epochs, epoch_t] = ecog_makeEpochs(data, newevents.event_sample, epochTime, srate);  
%[epochs] = ecog_normalizeEpochs(epochs, epoch_t, [-0.2 0], 'subtractwithintrial');

%% Compute spectra using Welch

fft_w    = window(@hann,500); % window for fft
fft_ov   = 100; % window overlap
reg_erp  = 0; 
t        = [0 0.5]; 
fft_t    = epoch_t > t(1) & epoch_t < t(2); 

% Dora's ecog_spectra function expects channels X epochs X time
data_epoch = permute(epochs, [3 2 1]);
stims = ones(1,size(data_epoch,2));

% Compute spectra
[f,spectra] = ecog_spectra(data_epoch,stims,fft_w,fft_t,fft_ov,srate,reg_erp);

%% Plot spectra

trial_ind = contains(events.task_name, 'prf') & ~contains(events.trial_name, 'BLANK');
blank_ind = contains(events.task_name, 'prf') & contains(events.trial_name, 'BLANK');

nChan = length(chan_inx); nSubPlot = ceil(sqrt(nChan));

f_ind = 2:200;

figure;
for ii = 1:nChan
    subplot(3, 2,ii); hold on
    chan_spectra = squeeze(spectra(chan_inx(ii),:,:));
    p1 = plot(f(f_ind), chan_spectra(trial_ind,f_ind), 'Color', [1 0.8 0.8], 'LineWidth', 1, 'LineStyle', ':');
    for jj = 1:length(p1), p1(jj).Annotation.LegendInformation.IconDisplayStyle = 'off'; end
    p2 = plot(f(f_ind), chan_spectra(blank_ind,f_ind), 'Color', [0.8 0.8 0.8], 'LineWidth', 1, 'LineStyle', ':');
    for jj = 1:length(p2), p2(jj).Annotation.LegendInformation.IconDisplayStyle = 'off'; end
    p1_m = plot(f(f_ind), mean(chan_spectra(trial_ind,f_ind)), 'r', 'LineWidth', 3);
    p2_m = plot(f(f_ind), mean(chan_spectra(blank_ind,f_ind)), 'k:', 'LineWidth', 3);
    if ii == 1, legend({'BAR', 'BLANK'});end
    set(gca, 'YScale', 'log', 'YLim', [10^-2.5 10^2.5]);
    title(channels.name(chan_inx(ii)));
    xlabel('frequency (Hz)'); ylabel('power spectral density estimate');
end
set(gcf, 'Position', [1000 500 1000 1000])

%% generate PRF time series

%f_ind = [30:49 51:99 101:149]; % skip 100 Hz
%f_ind = [30:45 55:95 105:145 155:200]; % skip 100 Hz
f_ind = [30:49 51:99 101:149 151:200]; % skip 100 Hz
%f_ind = [30:200];
prf_ts = geomean(spectra(:, :, f_ind),3); % take geomean to prevent bias to lower frequencies
nEvents = height(events);
prf_ts = reshape(prf_ts, [size(prf_ts,1) nEvents/2 2]);

figure;
for ii = 1:nChan
    subplot(3,2,ii); hold on
    plot(prf_ts(chan_inx(ii),:,1), 'r', 'LineWidth', 1)
    plot(prf_ts(chan_inx(ii),:,2), 'b', 'LineWidth', 1)
    plot(mean(prf_ts(chan_inx(ii),:,:),3), 'k:', 'LineWidth', 3)
    if ii == 1
        legend({'run1', 'run2', 'mean (run1 run2)'});
    end
    axis tight
    title(channels.name(chan_inx(ii)));
    xlabel('PRF stimulus index (bar position)'); ylabel(sprintf('average power between %s - %s Hz', num2str(min(f_ind)), num2str(max(f_ind))));
end
set(gcf, 'Position', [1000 500 1000 1000])

%% fit PRF model

bar_apertures = double(imresize(a.bar_apertures, [100 100], 'nearest'));

% Inputs to analyzePRF
stimulus = {bar_apertures};
data = {mean(prf_ts(chan_inx,:,:),3)}; % average over repeats
tr = 1;

opt = [];
opt.hrf = 1;
opt.maxpolydeg = 0;
opt.xvalmode = 0; 
opt.display = 'off';
opt.typicalgain = 1;
%opt.seedmode = -2;

% Run analyzePRF
results = analyzePRF_bounds(stimulus,data,tr,opt); 
chanNames = channels.name(chan_inx);

% Run analyzePRF bootstrap
nboots = 100;

clear results_boot 
for ii = 1:nboots 
    idx = randi(length(data{1}), [1 size(data{1},2)]);
    if ii == 1
        results_boot = analyzePRF_bounds(stimulus{1}(:,:,idx),data{1}(:,idx),tr,opt); 
    else 
        results_boot(ii) = analyzePRF_bounds(stimulus{1}(:,:,idx),data{1}(:,idx),tr,opt); 
    end
end

%% Save the results
save(fullfile(savePth, 'chaam_analyzePRF_results_ecog'), 'results', 'results_boot', 'chanNames');



%% %%%%%%%%%%%%%%%%%%


%% Visualize the location of each voxel's pRF

% The stimulus is 100 pixels (in both height and weight), and this corresponds to
% 16.6 degrees of visual angle:
cfactor = 16.6/100;

% Plot each electrode in separate subplots:
figure; hold on;

for cc = 1:nChan
        
    subplot(3,2,cc); hold on
    set(gcf,'Units','points','Position',[100 100 400 400]);

    xpos = results.ecc(cc) * cos(results.ang(cc)/180*pi) * cfactor;
    ypos = results.ecc(cc) * sin(results.ang(cc)/180*pi) * cfactor;
    ang = results.ang(cc)/180*pi;
    sd = results.rfsize(cc) * cfactor;
    h = k_drawellipse(xpos,ypos,ang,sd,sd);  % 
    set(h,'Color','k','LineWidth',2);
    set(scatter(xpos,ypos,'r.'),'CData',[0 0 0]);
    rsq = results.R2(cc);

    % plot edits
    h = k_drawellipse(0,0,0,8.3,8.3); % circle indicating stimulus extent
    set(h,'Color',[0.5 0.5 0.5],'LineWidth',1, 'LineStyle', ':');
    axis([-20 20 -20 20]);
    straightline(0,'h','k:');       % line indicating horizontal meridian
    straightline(0,'v','k:');       % line indicating vertical meridian
    axis square;
    set(gca,'XTick',-20:2:20,'YTick',-20:2:20);
    xlabel('X-position (deg)');
    ylabel('Y-position (deg)');

    title(sprintf('%s R2 = %s', chanNames{cc}, num2str(round(rsq,2))));
end
set(gcf, 'Position', [1000 500 1000 1000])

% All in one plot 

figure; hold on;
set(gcf,'Units','points','Position',[100 100 400 400]);
cmap = hsv(nChan);
for cc = 1:nChan 
      xpos = results.ecc(cc) * cos(results.ang(cc)/180*pi) * cfactor;
      ypos = results.ecc(cc) * sin(results.ang(cc)/180*pi) * cfactor;
      ang = results.ang(cc)/180*pi;
      sd = results.rfsize(cc) * cfactor;
      h = k_drawellipse(xpos,ypos,ang,sd,sd);  % 
      h.Annotation.LegendInformation.IconDisplayStyle = 'off';
      set(h,'Color',cmap(cc,:),'LineWidth',2);
      set(scatter(xpos,ypos,'r.'),'CData',cmap(cc,:));
end
h = k_drawellipse(0,0,0,8.3,8.3); % circle indicating stimulus extent
set(h,'Color',[0.5 0.5 0.5],'LineWidth',1, 'LineStyle', ':');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
axis([-20 20 -20 20]);
straightline(0,'h','k-');       % line indicating horizontal meridian
straightline(0,'v','k-');       % line indicating vertical meridian
axis square;
set(gca,'XTick',-20:2:20,'YTick',-20:2:20);
xlabel('X-position (deg)');
ylabel('Y-position (deg)');
legend(channels.name(chan_inx));
set(gca, 'FontSize', 18);
set(gcf, 'Position', [163   553   684   599]);

% Alternative visualization: plot the gaussian itself
% [xx, yy] = meshgrid(linspace(-1,1,100));
% [th, r] = cart2pol(xx, yy);
% mask = r < 1;
% 
% prms = NaN(20,5);
% imall = NaN(100,100,20);
% for ii =1:numel(r1)
%     p = r1(ii).params(1,:,63);
%     prms(ii,:) = p;
%     im = makegaussian2d(100,p(1),p(2),p(3)/sqrt(p(5)),p(3)/sqrt(p(5)));
%     im = im .* mask;
%     im = im / sum(im(:));
%     imall(:,:,ii) = im;
%      % makegaussian2d(res,r,c,sr,sc,xx,yy,ang,omitexp)
% end

%% Visualize the location of each voxel's pRF BOOTSTRAPS

% The stimulus is 100 pixels (in both height and weight), and this corresponds to
% 16.6 degrees of visual angle:
cfactor = 16.6/100;
clear xpos ypos ang sd rsq

% Plot each electrode in separate subplots:
figure; hold on;

for cc = 1:nChan
        
    subplot(3,2,cc); hold on
    set(gcf,'Units','points','Position',[100 100 400 400]);

    for ii = 1:nboots 

        theseresults = results_boot(ii); 
        
        xpos(cc,ii) = theseresults.ecc(cc) * cos(theseresults.ang(cc)/180*pi) * cfactor;
        ypos(cc,ii)= theseresults.ecc(cc) * sin(theseresults.ang(cc)/180*pi) * cfactor;
        ang(cc,ii) = theseresults.ang(cc)/180*pi;
        sd(cc,ii) = theseresults.rfsize(cc) * cfactor;
        h = k_drawellipse(xpos(cc,ii),ypos(cc,ii),ang(cc,ii),sd(cc,ii),sd(cc,ii));  
        set(h,'Color',[0.5 0.5 0.5],'LineStyle', '-','LineWidth',1);
        scatter(xpos(cc,ii),ypos(cc,ii),10,[0.5 0.5 0.5], 'o', 'filled');
        rsq(cc,ii) = theseresults.R2(cc);
    end
    
    title(sprintf('%s median R2 = %s', channels.name{chan_inx(cc)}, num2str(round(median(rsq(cc,:),2),1))));
    h = k_drawellipse(median(xpos(cc,:),2),median(ypos(cc,:),2),median(ang(cc,:),2),median(sd(cc,:),2),median(sd(cc,:),2)); 
    set(h,'Color','k','LineStyle', '-','LineWidth',3);
    h = scatter(median(xpos(cc,:),2),median(ypos(cc,:),2),'ko', 'filled');
    set(h, 'MarkerEdgeColor', 'k', 'SizeData',10);
    
    % plot edits
    h = k_drawellipse(0,0,0,8.3,8.3); % circle indicating stimulus extent
    set(h,'Color',[0.5 0.5 0.5],'LineWidth',1, 'LineStyle', ':');
    axis([-20 20 -20 20]);
    straightline(0,'h','k:');       % line indicating horizontal meridian
    straightline(0,'v','k:');       % line indicating vertical meridian
    axis square;
    set(gca,'XTick',-20:2:20,'YTick',-20:2:20);
    xlabel('X-position (deg)');
    ylabel('Y-position (deg)');

end
set(gcf, 'Position', [1000 500 1000 1000])

%% All in one plot  BOOTSTRAPS

figure; hold on;
set(gcf,'Units','points','Position',[100 100 400 400]);
cmap = hsv(nChan);
for cc = 1:nChan 
      h = k_drawellipse(median(xpos(cc,:),2),median(ypos(cc,:),2),median(ang(cc,:),2),median(sd(cc,:),2),median(sd(cc,:),2));
      h.Annotation.LegendInformation.IconDisplayStyle = 'off';
      set(h,'Color',cmap(cc,:),'LineWidth',2);
      set(scatter(median(xpos(cc,:),2),median(ypos(cc,:),2),'r.'),'CData',cmap(cc,:));
      set(h,'Color',cmap(cc,:));
end
h = k_drawellipse(0,0,0,8.3,8.3); % circle indicating stimulus extent
set(h,'Color',[0.5 0.5 0.5],'LineWidth',1, 'LineStyle', ':');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
axis([-20 20 -20 20]);
straightline(0,'h','k-');       % line indicating horizontal meridian
straightline(0,'v','k-');       % line indicating vertical meridian
axis square;
set(gca,'XTick',-20:2:20,'YTick',-20:2:20);
xlabel('X-position (deg)');
ylabel('Y-position (deg)');
legend(channels.name(chan_inx));
set(gca, 'FontSize', 18);
set(gcf, 'Position', [163   553   684   599]);

% %% Plot R2
% figure;plot(results.R2, 'LineWidth', 2);
% set(gca, 'FontSize', 18, 'XTick', 1:8:height(channels));
% xlabel('electrode inx'); ylabel('R2'); ylim([0 100]);
% 
% figure;histogram(results.R2, 100, 'FaceColor', 'b');
% set(gca, 'XLim', [0 100], 'FontSize', 18);
% xlabel('R2'); ylabel('number of elecs'); 

%% KK example code: plot time courses with fit

% Define some variables
res = [100 100];                    % row x column resolution of the stimuli
resmx = 100;                        % maximum resolution (along any dimension)
hrf = results.options.hrf;          % HRF that was used in the model
degs = results.options.maxpolydeg;  % vector of maximum polynomial degrees used in the model

% Pre-compute cache for faster execution
[d,xx,yy] = makegaussian2d(resmx,2,2,2,2);

% Prepare the stimuli for use in the model
stimulusPP = {};
for cc=1:length(stimulus)
  stimulusPP{cc} = squish(stimulus{cc},2)';  % this flattens the image so that the dimensionality is now frames x pixels
  stimulusPP{cc} = [stimulusPP{cc} cc*ones(size(stimulusPP{cc},1),1)];  % this adds a dummy column to indicate run breaks
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
for cc=1:length(degs)
  polymatrix{cc} = projectionmatrix(constructpolynomialmatrix(size(data{cc},2),0:degs(cc)));
end

% Visualize the results
figure; hold on;

for ii=1:nChan 
    
    vx = chan_inx(ii);
    % For each run, collect the data and the model fit.  We project out polynomials
    % from both the data and the model fit.  This deals with the problem of
    % slow trends in the data.
    datats = {};
    modelts = {};
    for cc=1:length(data)
      datats{cc} =  polymatrix{cc}*data{cc}(vx,:)';
      modelts{cc} = polymatrix{cc}*modelfun(results.params(1,:,vx),stimulusPP{cc});
    end
 
    subplot(3, 2,ii); hold on
    set(gcf,'Units','points');
    plot(cat(1,datats{:}),'k-', 'LineWidth', 2);
    plot(cat(1,modelts{:}),'r-','LineWidth', 2);
    xlabel('PRF stimulus','FontSize', 28);
    ax = axis;
    %axis([.5 1200+.5 ax(3:4)]);
    if ii == 1, legend('Data', 'Model prediction'); end
    set(gca, 'FontSize', 18)
    title(channels.name(chan_inx(ii)));
    axis tight
end

%% Comparison with benson

% The stimulus is 100 pixels (in both height and weight), and this corresponds to
% 16 degrees of visual angle.  To convert from pixels to degreees, we multiply
% by 10/100.
cfactor =16.6/100;
prf_ecc = results.ecc*cfactor;
thresh = 30;

figure;hold on;

ind = results.R2>thresh; %& results.ecc < 20;
r = corr(prf_ecc(ind), channels.bensoneccen(ind),'rows','complete', 'Type', 'Spearman');
subplot(2,1,1);scatter(channels.bensoneccen(ind),prf_ecc(ind), 100, 'filled');
ylabel('eccentricity from analyzePRF');
xlabel('eccentricity from benson template');
set(gca, 'FontSize', 18); 
title(sprintf('for R2 > %d, correlation = %s',thresh, num2str(round(r,2)))); 

%ang = ang + 90;%ang = mod(90-(-ang),360);
ind = results.R2>thresh; %&results.ecc < 50;
r = corr(results.ang(ind), channels.bensonangle(ind)+90,'rows','complete',  'Type', 'Spearman');
subplot(2,1,2);scatter(channels.bensonangle(ind)+90,results.ang(ind),100,'filled');
ylabel('angle from analyzePRF');
xlabel('angle from benson template');
axis([0 360 0 360]);
set(gca, 'FontSize', 18);
title(sprintf('for R2 > %d, correlation = %s', thresh, num2str(round(r,2)))); 

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
set(gca, 'FontSize', 18);

