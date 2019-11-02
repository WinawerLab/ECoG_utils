
%% Paths specs

% Input path
dataPth     = '/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/ECoGCAR';

% Dataset specs
subject     = 'chaam'; 
session     = 'umcuecog';
task        = {'prf'};
description = 'reref';

% Get bar apertures (necessary for PRF fitting)
%a = load('/Volumes/GoogleDrive/My Drive/temp/bar_apertures.mat');
a = load('/Users/winawerlab/matlab/toolboxes/BAIRstimuli/stimuli/bar_apertures.mat');

% Channels of interest
COI = {'Oc17','Oc18', 'Oc24', 'Oc16', 'sT1', 'Oc31'}; 

%% Get the data
fprintf('Reading in the data...\n');
[data, channels, events] = bidsEcogGetPreprocData(dataPth, subject, session, task, [], description);
srate = channels.sampling_frequency(1);
chan_inx = find(contains(channels.name, COI));

%% Get visual area matches for this subject

fprintf('Computing matches with visual atlases...\n');
specs = [];
specs.pID           = subject; 
specs.plotmesh      = 'right';
specs.plotlabel     = 'yes'; %specs.atlasNames    = {'benson14_varea'};
BIDSformatted       = 1;
visualelectrodes    = electrode_to_nearest_node(specs, BIDSformatted);

% Add visual area names (W and B) ecc, angle, sigma to channels table
[channels] = bair_addVisualAtlasNamesToChannelTable(channels,visualelectrodes);

%% Data preprocessing

% % Shift onsets (estimated delay with respect to NYU data)
fprintf('[%s] This is a umcu patient: shifting the onsets\n',mfilename);
shiftInSeconds = 0.062; % 62 ms
shiftInSamples = round(shiftInSeconds*srate); 
events.onset = events.onset + shiftInSeconds;
events.event_sample = events.event_sample + shiftInSamples; 

% Epoch the data
epochTime = [-0.3 0.85];
fprintf('Computing epochs...\n');
[epochs, epoch_t] = ecog_makeEpochs(data, events.event_sample, epochTime, srate);  
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

figure;
for ii = 1:nChan
    subplot(3, 2,ii); hold on
    chan_spectra = squeeze(spectra(chan_inx(ii),:,:));
    plot(f(f_ind), mean(chan_spectra(trial_ind,f_ind)) - mean(chan_spectra(blank_ind,f_ind)), 'k', 'LineWidth', 2);
    if ii == 1, legend({'BAR-BLANK'}); end
    l1 = line([f(f_ind(1)) f(f_ind(end))],[0 0], 'Color', 'r', 'LineStyle', ':', 'LineWidth', 2);
    l1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    set(gca, 'YLim', [-30 30]);
    title(channels.name(chan_inx(ii)));
    xlabel('frequency (Hz)'); ylabel('power spectral density estimate difference');

end

%% generate PRF time series

f_ind = [30:49 51:99 101:149]; % skip 100 Hz
prf_ts = geomean(spectra(:, :, f_ind),3); % take geomean to prevent bias to lower frequencies
nEvents = height(events);
prf_ts = reshape(prf_ts, [size(prf_ts,1) nEvents/2 2]);

figure;
for ii = 1:nChan
    subplot(3,2,ii); hold on
    plot(prf_ts(chan_inx(ii),:,1), 'r', 'LineWidth', 1)
    plot(prf_ts(chan_inx(ii),:,2), 'b', 'LineWidth', 1)
    plot(mean(prf_ts(chan_inx(ii),:,:),3), 'k', 'LineWidth', 2)
    if ii == 1
        legend({'run1', 'run2', 'mean (run1 run2)'});
    end
    axis tight
    title(channels.name(chan_inx(ii)));
    xlabel('PRF stimulus index (bar position)'); ylabel(sprintf('average power between %s - %s Hz', num2str(min(f_ind)), num2str(max(f_ind))));
end

%% fit PRF model

bar_apertures = imresize(a.bar_apertures, [100 100], 'nearest');

% Inputs to analyzePRF
stimulus = {bar_apertures};
data = {mean(prf_ts,3)}; % average over repeats
tr = 1;

opt.hrf = 1;
opt.maxpolydeg = 0;
opt.xvalmode = 0; 
opt.display = 'off';

% Run analyzePRF
results = analyzePRF(stimulus,data,tr,opt); 

% % Run analyzePRF bootstrap
% nboots = 100;
% clear r1 
% for ii = 1:nboots 
%     idx = randi(length(data{1}), [1 size(data{1},2)]);
%     if ii == 1
%         r1 = analyzePRF(stimulus{1}(:,:,idx),data{1}(:,idx),tr,opt); 
%     else 
%         r1(ii) = analyzePRF(stimulus{1}(:,:,idx),data{1}(:,idx),tr,opt); % WITHOUT baseline correction
%     end
% end

%% Plot R2
figure;plot(results.R2, 'LineWidth', 2);
set(gca, 'FontSize', 18, 'XTick', 1:8:height(channels));
xlabel('electrode inx'); ylabel('R2'); ylim([0 100]);

figure;histogram(results.R2, 100, 'FaceColor', 'b');
set(gca, 'XLim', [0 100], 'FontSize', 18);
xlabel('R2'); ylabel('number of elecs'); 

%% %% KK example code: Visualize the location of each voxel's pRF
% The stimulus is 100 pixels (in both height and weight), and this corresponds to
% 16.6 degrees of visual angle.  To convert from pixels to degreees, we multiply
% by 16.6/100.
cfactor = 16.6/100;
clear xpos ypos ang sd rsq

figure; hold on;

for cc = 1:nChan
        
    subplot(nSubPlot, nSubPlot,cc); hold on
    set(gcf,'Units','points','Position',[100 100 400 400]);

    for ii = 1:nboots 

        results = r1(ii); 
        xpos(cc,ii) = results.ecc(chan_inx(cc)) * cos(results.ang(chan_inx(cc))/180*pi) * cfactor;
        ypos(cc,ii) = results.ecc(chan_inx(cc)) * sin(results.ang(chan_inx(cc))/180*pi) * cfactor;
        ang(cc,ii) = results.ang(chan_inx(cc))/180*pi;
        sd(cc,ii) = results.rfsize(chan_inx(cc)) * cfactor;
        h = k_drawellipse(xpos(cc,ii),ypos(cc,ii),ang(cc,ii),2*sd(cc,ii),2*sd(cc,ii));  
        set(h,'Color',[0.5 0.5 0.5],'LineStyle', '-','LineWidth',1);
        scatter(xpos(cc,ii),ypos(cc,ii),10,[0.5 0.5 0.5], 'o', 'filled');
        rsq(cc,ii) = results.R2(chan_inx(cc));

        % plot edits
        k_drawrectangle(0,0,16.6,16.6,'k-');  % square indicating stimulus extent
        axis([-20 20 -20 20]);
        straightline(0,'h','k-');       % line indicating horizontal meridian
        straightline(0,'v','k-');       % line indicating vertical meridian
        axis square;
        set(gca,'XTick',-20:2:20,'YTick',-20:2:20);
        xlabel('X-position (deg)');
        ylabel('Y-position (deg)');
    end

    title(sprintf('%s median R2 = %s', channels.name{chan_inx(cc)}, num2str(round(median(rsq(cc,:),2),1))));
    h = k_drawellipse(median(xpos(cc,:),2),median(ypos(cc,:),2),median(ang(cc,:),2),median(2*sd(cc,:),2),median(2*sd(cc,:),2)); 
    set(h,'Color','k','LineStyle', '-','LineWidth',3);
    h = scatter(median(xpos(cc,:),2),median(ypos(cc,:),2),'ko', 'filled');
    set(h, 'MarkerEdgeColor', 'k', 'SizeData',10);
end

%% %% KK example code: Visualize the location of each voxel's pRF

% The stimulus is 100 pixels (in both height and weight), and this corresponds to
% 16.6 degrees of visual angle.  To convert from pixels to degreees, we multiply
% by 16.6/100.
cfactor = 16.6/100;

figure; hold on;
set(gcf,'Units','points','Position',[100 100 400 400]);
cmap = jet(nChan);
for p=1:nChan %size(results.ang,1)
      xpos = results.ecc(chan_inx(p)) * cos(results.ang(chan_inx(p))/180*pi) * cfactor;
      ypos = results.ecc(chan_inx(p)) * sin(results.ang(chan_inx(p))/180*pi) * cfactor;
      ang = results.ang(chan_inx(p))/180*pi;
      sd = results.rfsize(chan_inx(p)) * cfactor;
      h = k_drawellipse(xpos,ypos,ang,sd,sd);  % 
      h.Annotation.LegendInformation.IconDisplayStyle = 'off';
      set(h,'Color',cmap(p,:),'LineWidth',2);
      set(scatter(xpos,ypos,'r.'),'CData',cmap(p,:));
end
k_drawrectangle(0,0,16.6,16.6,'k-');  % square indicating stimulus extent
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

% Visualize the results
figure; hold on;

for ii=1:nChan 
    
    vx = chan_inx(ii);
    % For each run, collect the data and the model fit.  We project out polynomials
    % from both the data and the model fit.  This deals with the problem of
    % slow trends in the data.
    datats = {};
    modelts = {};
    for p=1:length(data)
      datats{p} =  polymatrix{p}*data{p}(vx,:)';
      modelts{p} = polymatrix{p}*modelfun(results.params(1,:,vx),stimulusPP{p});
    end
 
    subplot(nSubPlot, nSubPlot,ii); hold on
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

