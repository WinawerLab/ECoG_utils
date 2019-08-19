
tbUse ECoG_utils;

%% Define paths and dataset

% Dataset specs
projectName = 'visual';
sub_label   = 'som726'; 
ses_label   = {'nyuecog02', 'nyuecog03'};

% Input paths 
projectDir  = '/Volumes/server/Projects/BAIR/';
dataPth     = fullfile(projectDir, 'Data', 'BIDS');

%% Load preprocessed data

% Output paths specs
for ii = 1:length(ses_label)
    procDir = fullfile(dataPth, projectName, 'derivatives', 'preprocessed', sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label{ii}));
    dataName = fullfile(procDir, sprintf('sub-%s_ses-%s_epoched.mat', sub_label, ses_label{ii}));
    all{ii} = load(dataName);
end

% Concatenate the sessions
trials = all{1}.trials;
trials.events = [all{1}.trials.events; all{2}.trials.events];
trials.broadband = cat(3,all{1}.trials.broadband, all{2}.trials.broadband);

%% Select data from a single channel

% Select a subset of electrodes to analyze
elecIndex  = contains(trials.channels.name, 'GB119'); % electrode with good signal
channel   = trials.channels(elecIndex,:);

% Select the prf events
trialIndex = contains(trials.events.task_name, 'prf');
events     = trials.events(trialIndex,:);

% Extract the prf data
bb         = trials.broadband(elecIndex,:,trialIndex);

%% Perform baseline correction

% % From ZHOU et al 2019: To convert the unit of the time-varying broadband
% to percent signal change in each electrode, we first averaged each
% broadband time series across epochs. We defined the first 200 ms prior to
% stimulus onset as the baseline period for the epoch-averaged time course,
% then we computed the percent signal change by dividing the entire 1200ms
% time course point-wise by the average of the baseline. To equalize the
% baseline across electrodes, we subtracted the baseline average from the
% entire time-course so each electrode has trial-averaged baseline 0.
base_range = (trials.time >= -0.2 & trials.time < 0);
m_base     = squeeze(median(mean(bb(:, base_range, :), 2), 3));
bbWBC      = bb./m_base-1;

%% Compute response magnitude per stimulus

% Define a time window over which to average the broadband timecourse
time_win   = [0.05 0.55];

% Compute average broadband response in time window

% WITHOUT baseline correction:
ts1       = squeeze(mean(bb(:,trials.time>time_win(1) & trials.time<time_win(2),:),2));  
% WITH baseline correction:
ts2       = squeeze(mean(bbWBC(:,trials.time>time_win(1) & trials.time<time_win(2),:),2)); 

% Plot the PRF response
figure;hold on;
plot(ts1(:,:,1))
plot(ts2(:,:,1))
legend('without baseline correction', 'with baseline correction');
corr(ts1,ts2)

%% Perform PRF fits

% Load apertures
load('/Users/winawerlab/matlab/toolboxes/BAIRstimuli/stimuli/bar_apertures.mat');
bar_apertures = imresize(bar_apertures, [100 100], 'nearest');

% Inputs to analyzePRF
stimulus = {bar_apertures,bar_apertures,bar_apertures,bar_apertures}; % 4 runs
tr = 1;
opt.hrf = 1;
opt.maxpolydeg = 0; % polynomial of degree 0
opt.xvalmode = 0; % fit all the data

% Separate the four runs: reshape then put into cell (input for analyzePRF)
data1 = reshape(ts1,[size(bb,1) size(bb,3)/4 4]); % without baseline correction
data2 = reshape(ts2,[size(bb,1) size(bb,3)/4 4]); % with baseline correction
data1 = {data1(:,:,1),data1(:,:,2),data1(:,:,3),data1(:,:,4)}; % without baseline correction
data2 = {data2(:,:,1),data2(:,:,2),data2(:,:,3),data2(:,:,4)}; % with baseline correction

% Run analyzePRF
results1 = analyzePRF(stimulus,data1,tr,opt); % WITHOUT baseline correction
results2 = analyzePRF(stimulus,data2,tr,opt); % WITH baseline correction

%% KK example code: plot time courses with fit

% Define some variables
res = [100 100];                    % row x column resolution of the stimuli
resmx = 100;                        % maximum resolution (along any dimension)
hrf = results.options.hrf;          % HRF that was used in the model
degs = results1.options.maxpolydeg;  % vector of maximum polynomial degrees used in the model

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

% For each run, collect the data and the model fit.  We project out polynomials
% from both the data and the model fit.  This deals with the problem of
% slow trends in the data.
for p=1:length(data)
  datats1{p} =  polymatrix{p}*data1{p}(vx,:)';
  modelts1{p} = polymatrix{p}*modelfun(results1.params(1,:,vx),stimulusPP{p});
  datats2{p} =  polymatrix{p}*data2{p}(vx,:)';
  modelts2{p} = polymatrix{p}*modelfun(results2.params(1,:,vx),stimulusPP{p});
end

% Compare the data time series with and wo baseline correction
figure; hold on;
plot(cat(1,datats1{:}), 'k', 'LineWidth', 2)
plot(cat(1,datats2{:}), 'b', 'LineWidth', 2)
corr(cat(1,datats1{:}), cat(1,datats2{:}))
legend('without baseline correction', 'with baseline correction');
title('data');

% Compare the predicted time series with and wo baseline correction
figure; hold on;
plot(cat(1,modelts1{:}), 'k', 'LineWidth', 2)
plot(cat(1,modelts2{:}), 'b', 'LineWidth', 2)
corr(cat(1,modelts1{:}), cat(1,modelts2{:}))
legend('without baseline correction', 'with baseline correction');
title('model prediction');

% % Visualize the results
% figure; hold on;
% set(gcf,'Units','points','Position',[100 100 1000 100]);
% plot(cat(1,datats1{:}),'k-', 'LineWidth', 2);
% plot(cat(1,modelts1{:}),'r-','LineWidth', 2);
% straightline(224*(1:2)+.5,'v','g-');
% xlabel('PRF stimulus','FontSize', 28);
% ylabel('Broadband response','FontSize', 28);
% ax = axis;
% %axis([.5 1200+.5 ax(3:4)]);
% legend('data', 'model prediction');
% 
% set(gcf, 'Position', [60 300 2000 1000]);
% set(gca, 'FontSize', 18)
% 
% % Visualize the results
% figure; hold on;
% set(gcf,'Units','points','Position',[100 100 1000 100]);
% plot(cat(1,datats2{:}),'k-', 'LineWidth', 2);
% plot(cat(1,modelts2{:}),'r-','LineWidth', 2);
% straightline(224*(1:2)+.5,'v','g-');
% xlabel('PRF stimulus','FontSize', 28);
% ylabel('Broadband response','FontSize', 28);
% ax = axis;
% %axis([.5 1200+.5 ax(3:4)]);
% legend('data', 'model prediction');
% 
% set(gcf, 'Position', [60 300 2000 1000]);
% set(gca, 'FontSize', 18)

% %% %% KK example code: Visualize the location of each voxel's pRF
% 
% % The stimulus is 100 pixels (in both height and weight), and this corresponds to
% % 16.6 degrees of visual angle.  To convert from pixels to degreees, we multiply
% % by 16.6/100.
% cfactor = 16.6/100;
% 
% % WITHOUT baseline correction
% results = results1; 
% figure; hold on;
% set(gcf,'Units','points','Position',[100 100 400 400]);
% cmap = jet(size(results.ang,1));
% for p=1:size(results.ang,1)
%   if results.R2(p)>30
%       xpos = results.ecc(p) * cos(results.ang(p)/180*pi) * cfactor;
%       ypos = results.ecc(p) * sin(results.ang(p)/180*pi) * cfactor;
%       ang = results.ang(p)/180*pi;
%       sd = results.rfsize(p) * cfactor;
%       h = k_drawellipse(xpos,ypos,ang,2*sd,2*sd);  % circle at +/- 2 pRF sizes
%       set(h,'Color',cmap(p,:),'LineWidth',2);
%       set(scatter(xpos,ypos,'r.'),'CData',cmap(p,:));
%   end
% end
% k_drawrectangle(0,0,16.6,16.6,'k-');  % square indicating stimulus extent
% axis([-20 20 -20 20]);
% straightline(0,'h','k-');       % line indicating horizontal meridian
% straightline(0,'v','k-');       % line indicating vertical meridian
% axis square;
% set(gca,'XTick',-20:2:20,'YTick',-20:2:20);
% xlabel('X-position (deg)');
% ylabel('Y-position (deg)');
% 
% % WITH baseline correction
% results = results2; 
% figure; hold on;
% set(gcf,'Units','points','Position',[100 100 400 400]);
% cmap = jet(size(results.ang,1));
% for p=1:size(results.ang,1)
%   if results.R2(p)>30
%       xpos = results.ecc(p) * cos(results.ang(p)/180*pi) * cfactor;
%       ypos = results.ecc(p) * sin(results.ang(p)/180*pi) * cfactor;
%       ang = results.ang(p)/180*pi;
%       sd = results.rfsize(p) * cfactor;
%       h = k_drawellipse(xpos,ypos,ang,2*sd,2*sd);  % circle at +/- 2 pRF sizes
%       set(h,'Color',cmap(p,:),'LineWidth',2);
%       set(scatter(xpos,ypos,'r.'),'CData',cmap(p,:));
%   end
% end
% k_drawrectangle(0,0,16.6,16.6,'k-');  % square indicating stimulus extent
% axis([-20 20 -20 20]);
% straightline(0,'h','k-');       % line indicating horizontal meridian
% straightline(0,'v','k-');       % line indicating vertical meridian
% axis square;
% set(gca,'XTick',-20:2:20,'YTick',-20:2:20);
% xlabel('X-position (deg)');
% ylabel('Y-position (deg)');