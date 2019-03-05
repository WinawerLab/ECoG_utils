rootpath = '/Volumes/server/Projects/BAIR/Data/Raw/ECoG/';
sub_label = '661'; % channel with triggers = 138 (DC10_REF)

dataName = fullfile(rootpath, sub_label, ['NY' sub_label '_BairPilot.edf']);
stimFile = fullfile(rootpath, sub_label, ['NY' sub_label '_StimInfo.mat']);
%tbUse('fieldtrip');
%addpath(genpath([rootpath '/code']));

%% Import data

cfg            = [];
cfg.dataset    = dataName;
cfg.continuous = 'yes';
cfg.channel    = 'all';
data           = ft_preprocessing(cfg);

load(stimFile);

%% Remove parts of non-task data
% In this data set, there is a large section of recording prior to the
% start of the actual experiment:

% Plot one channel plus onsets:
figure;hold on
plot(data.time{1}, data.trial{1}(1,:)); 
plot(onsets.prf, zeros(length(onsets.prf),1), 'ok')
plot(onsets.soc, zeros(length(onsets.soc),1), 'ok')
title('All data');

% Removing some of these data to reduce file sizes:
removePre  = 1701; % seconds
removePost = 3800; % seconds
[~, startAt] = min(abs(data.time{1}-removePre));
[~, stopAt]  = min(abs(data.time{1}-removePost));

data.time{1}  = data.time{1}(startAt:stopAt);
data.trial{1} = data.trial{1}(:,startAt:stopAt);
data.sampleinfo = startAt:stopAt;
data.hdr.nSamples = length(startAt:stopAt);

figure;hold on
plot(data.time{1}, data.trial{1}(1,:)); 
plot(onsets.prf, zeros(length(onsets.prf),1), 'ok')
plot(onsets.soc, zeros(length(onsets.soc),1), 'ok')
title('Reduced data');

%% Look at powerpectrum of channels to see which ones look bad

% Look at data.label to determine which channels are EEG channels: 1:110
figure;spectopo(data.trial{1}(1:110,:),size(data.trial{1},2),data.fsample);

[pxx,f] = pwelch(data.trial{1}(1:110,:)',data.fsample,0,data.fsample,data.fsample);
figure,plot(f,log10(pxx))

% % % To inspect data on a channel by channel basis:
% for ii = 1:110
%     figure;hold on;
%     plot(data.time{1}, data.trial{1}(ii,:)); 
%     title([num2str(ii) ' ' data.label{ii}]);
%     % add onsets to indicate relevant data segments
%     plot(onsets.prf, zeros(length(onsets.prf),1), 'ok')
%     plot(onsets.soc, zeros(length(onsets.soc),1), 'ok')
%     waitforbuttonpress; close
% end
%     
bad_channels = [100 102 89]; % manually selected based on inspection of spectra and time courses
% 100 and 102 have partly flat spectra; 89 has line noise peak on 200 Hz

% Other observations:
% 33, 24, 48, 93 have a strong LF drift before ~stimulus onset first task
% 47, 57, 60 have a multiple of these LF drifts throughout tasks
% 69-71 have one (1) big spike during one of the PRF runs. Reason to remove?
% 105-110 show strong bursts of activity IN BETWEEN tasks!? Could be
% movement artifacts? Are these ECOG channels? --> checked with Hugh: these
% are 'subgaleal' electrodes that are not included in recon; exclude
good_channels = setdiff(1:104,bad_channels);

%% So let's visualize the timeseries of some channels
figure
subplot(2,1,1),hold on
plot(data.time{1},data.trial{1}(bad_channels(1),:),'r')
title('first bad channel')
subplot(2,1,2),hold on
plot(data.time{1},data.trial{1}(good_channels(1),:),'k')
title('first good channel')

% check the powerplot of all the good channels, no leftover outliers?
figure,plot(f,log10(pxx(:,good_channels)))

% check the timeseries of all the good channels, no noisy moments?
figure,plot(data.time{1},data.trial{1}(good_channels,:))

%%  Now we know the bad channels, and we can do a common average reference.
signal = data.trial{1};
[signal] = ecog_CarRegress(signal, good_channels);

%% look at the effect of CAR
figure
subplot(1,2,1),hold on
channel_plot = good_channels(1);
plot(data.time{1},data.trial{1}(channel_plot,:),'k')
plot(data.time{1},signal(channel_plot,:),'g')
legend({'before CAR','after CAR'})

subplot(1,2,2),hold on
[pxx2,f] = pwelch(signal',data.fsample,0,data.fsample,data.fsample);
plot(f,pxx(:,channel_plot),'k')
plot(f,pxx2(:,channel_plot),'g'); 
set(gca, 'YScale', 'log')

% Q what are the high frequency spectral peaks around 200 Hz? Do not all
% seem to be harmonics of line noise (60, 120, 180, 240)

%% Now you have preprocessed data, save it somewhere and analyse.
save(fullfile(rootpath, sub_label, sprintf('NY%s_data_preproc', sub_label)), 'data', 'signal', '-v7.3');

% NEXT: filter, then split data by task

%% You may want to do a notch filter (optional).

% Filter: Dora style; compare with other approaches?
% You can also notch filter at 100 and then take hband = 60:120

% hband_sig=zeros(size(signal));
% for el=1:size(signal,1)
%     disp(['filter el ' int2str(el) ' of ' int2str(size(signal,1))])
%     hband=[60 70];
%     hband_sig1=butterpass_eeglabdata(signal(el,:)',hband,data.fsample);   
%     hband_sig1=log10(abs(hilbert(hband_sig1)).^2)-mean(log10(abs(hilbert(hband_sig1)).^2));
%     hband=[70 80];
%     hband_sig2=butterpass_eeglabdata(signal(el,:)',hband,data.fsample);   
%     hband_sig2=log10(abs(hilbert(hband_sig2)).^2)-mean(log10(abs(hilbert(hband_sig2)).^2));
%     hband=[80 90];
%     hband_sig3=butterpass_eeglabdata(signal(el,:)',hband,data.fsample);   
%     hband_sig3=log10(abs(hilbert(hband_sig3)).^2)-mean(log10(abs(hilbert(hband_sig3)).^2));
%     hband=[110 120];
%     hband_sig4=butterpass_eeglabdata(signal(el,:)',hband,data.fsample);   
%     hband_sig4=log10(abs(hilbert(hband_sig4)).^2)-mean(log10(abs(hilbert(hband_sig4)).^2));
% 
%     hband_sig(el,:)=mean([hband_sig1 hband_sig2 hband_sig3 hband_sig4],2);
% 
%     clear hband_sig1 hband_sig2 hband_sig3 hband_sig4
% end

% compare with broadband computation function from Jons broadband tutorial
% bands = [[60 70]; [70 80]; [80 90]; [110 120]]; 
bands = [[70 90]; [90 110]; [130 150]; [150 170]];
% bands = [[70 80]; [80 90]; [90 100]; [100 110];[130 140]; [140 150]; [150 160]; [160 170]];
%                 1 'abs(hilbert(mean(whiten(bp(x)))))';
%                 2 'abs(hilbert(mean(whiten(bp(x))))).^2';
%                 3 'geomean(abs(hilbert(whiten(bp(x)))))';
%                 4 'geomean(abs(hilbert(whiten(bp(x))).^2)';
%                 5 'geomean(abs(hilbert(bp(x)).^2)';

broadband = ecog_extractBroadband(signal', data.fsample, 7, bands);
% BANDPASS FILTER THE SIGNAL
% bp = bandpassFilter(signal', data.fsample, bands);
% 
% % EXTRACT BROADBAND
% whiten         = @(x) (x - mean(x(:)))./ diff(prctile(x, [.25 .75]));
% whiten_hilbert = @(x) abs(hilbert(whiten(x)));
% 
% broadband = sum(whiten_hilbert(bp).^2, 3);
% 
hband_sig = broadband;

%% plot filtered time course for channel(s) of interest

% electrodes with coverage in visual areas: 
eltomatch = {'EEG OC_01-REF'}; % V3d (V3)
%eltomatch = {'OC_02'}; % V3b (V3)
%eltomatch = {'OC_05'}; % MT
%eltomatch = {'DPMT_01'}; % depth electrode matching w V2

% find matching electrode numbers
el = ecog_matchChannels(eltomatch, data);

load(stimFile);
tasklabels = fields(onsets);

for tl = 1:length(tasklabels)
    tasklabel = tasklabels{tl};
    
    % smooth & plot
    figure('Name', tasklabel), hold on,
    for ii = 1:length(eltomatch)
        plot(data.time{1},smooth(hband_sig(el(ii),:),128), 'LineWidth', 2);
        %plot(data.time{1},smooth(signal(el(ii),:),128), 'LineWidth', 2);
    end

    % and plot event onsets on top
    plot(onsets.(tasklabel), zeros(length(onsets.(tasklabel)),1),'k.','MarkerSize', 25, 'LineStyle','none');
    legend(data.label(el));
    xlabel('time (s)');
    ylabel('broadband power (70-170 Hz)');
    %ylabel('evoked');
    l = line([data.time{1}(1) max(data.time{1})], [0 0],'LineStyle', '-', 'Color', 'k');
    l.Annotation.LegendInformation.IconDisplayStyle = 'off';
end

%% epoch data (all channels, all tasks)

% set epoch length (in seconds)
prestim = 0.5; 
poststim = 1.5;

% determine how many samples to go back and forth to extract epoch
onset_pre = round(prestim/(1/data.fsample));
onset_post = round(poststim/(1/data.fsample));

trials = struct();
for tl = 1:length(tasklabels)
    
    tasklabel = tasklabels{tl};
    [~,onsetsInx] = intersect(data.time{1},onsets.(tasklabel));
    % for prf, skip second onsets
    if strmatch(tasklabel, 'prf')
        onsetsInx = onsetsInx(1:2:end);
    end

    % extract epochs
    for ii = 1:length(data.label)
        for jj = 1:length(onsetsInx)  
            % broadband
            trials.broadband.(tasklabel)(ii,:,jj) = smooth(hband_sig(ii,onsetsInx(jj)-onset_pre:onsetsInx(jj)+onset_post),8);
            % evoked
            trials.evoked.(tasklabel)(ii,:,jj) = smooth(signal(ii,onsetsInx(jj)-onset_pre:onsetsInx(jj)+onset_post),8);
        end
    end
end
disp('done');

trials.label = data.label;
trials.time = -prestim:(1/data.fsample):poststim;
trials.stimuli = stimuli;
trials.fsample = data.fsample;

%% Now you have epoched data, save it somewhere and analyse.
savePth = [rootpath sub_label];
save(fullfile(rootpath, sub_label, sprintf('NY%s_data_epoched_bbmethod7', sub_label)), 'trials', 'category_names', '-v7.3');

%% plot trial averages for channel(s) of interest

usebootstrap = 'no'; % yes/no
signaltoplot = 'broadband'; % evoked/broadband

% electrodes with coverage in visual areas: 
eltomatch = {'EEG OC_01-REF','EEG OC_02-REF', 'EEG OC_03-REF', 'EEG OC_04-REF', 'EEG OC_05-REF', 'EEG OC_05-REF'};

% find matching electrode numbers
el = ecog_matchChannels(eltomatch, trials);

for tl = 1:length(tasklabels)
    
    tasklabel = tasklabels{tl};
    ga = nan(length(el), length(trials.time)); 
    ulim = nan(length(el), length(trials.time));  
    llim = nan(length(el), length(trials.time)); 
    
    switch usebootstrap
        case 'yes'       
            % bootstrap to get distribution median and 95% intervals
            for ee = 1:length(el)
                for trial = 1:length(trials.time)
                    [bootstat] = bootstrp(1000,@mean,squeeze(trials.(signaltoplot).(tasklabel)(el(ee),trial,:)));
                    ga(ee,trial) = median(bootstat);
                    llim(ee,trial) = quantile(sort(bootstat), 0.05);
                    ulim(ee,trial) = quantile(sort(bootstat), 0.95);
                end
            end          
        otherwise
            % take average and standard deviation across trials
            ga = mean(trials.(signaltoplot).(tasklabel)(el,:,:),3);
            llim = ga-std(trials.(signaltoplot).(tasklabel)(el,:,:),0,3);
            ulim = ga+std(trials.(signaltoplot).(tasklabel)(el,:,:),0,3);
    end
    
    % plot   
    figure('Name', tasklabel); hold on;
    colors = jet(length(el));
       
    for ee = 1:length(el)
        h = ciplot(llim(ee,:),ulim(ee,:),trials.time,colors(ee,:),0.1);
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        plot(trials.time, ga(ee,:),'Color',colors(ee,:), 'LineWidth',2);
    end
    
    % determine axis limits
    lim = max(max(abs(ga)));
    ylim = [-lim-0.2*lim lim+0.2*lim];
    set(gca, 'YLim', ylim);
    
    % add stim onset and zero lines
    line([0 0], ylim,'LineStyle', ':', 'Color', 'k');
    line([-0.5 1.5], [0 0],'LineStyle', ':', 'Color', 'k');
    
    % add legends, labels
    legend(data.label(el));
    xlabel('time (s)');
    title(tasklabel);
    axis tight
    ylabel(signaltoplot);
end

%% plot individual trials: prf

tasklabel = 'prf'; 
signaltoplot = 'broadband'; 
eltoplot = 'EEG OC_01-REF';
%eltoplot = {'EEG OC_01-REF','EEG OC_02-REF', 'EEG OC_03-REF', 'EEG OC_04-REF', 'EEG OC_05-REF', 'EEG OC_05-REF'};

el = ecog_matchChannels(eltoplot,trials);
temp = squeeze(trials.(signaltoplot).(tasklabel)(el,:,:)); 

% plot trials in separate panels
figure('Name', [tasklabel ' ' signaltoplot ' ' trials.label{el}]); hold on
lim = max(max(abs(temp)));
ylim = [-lim-0.1*lim lim+0.1*lim];
ntrial = size(temp,2)/4;

% for trial = 1:ntrial
%     subplot(ceil(sqrt(ntrial)),ceil(sqrt(ntrial)),trial); hold on;
%     
%     plot(trials.time, temp(:,trial), 'k', 'LineWidth',2);
%     set(gca, 'YLim', ylim);
%     line([0 0], ylim,'LineStyle', ':', 'Color', 'k');
%     line([-0.1 0.5], [0 0],'LineStyle', ':', 'Color', 'k');
%     xlabel('time (s)');
%     %ylabel(ylabels{signaltoplot})
%     axis tight
%     if trial == 1
%         title([{'trial number:'} {num2str(trial)}])
%     else
%         title(num2str(trial));
%     end
% end

% plot trials on top of one another
figure('Name', [tasklabel ' ' signaltoplot ' ' trials.label{el}]); hold on
plot(trials.time,temp);
switch signaltoplot
    case 'broadband'
        ylim = [0 lim+0.2*lim];
    case 'evoked'
        ylim = [-lim-0.2*lim lim+0.2*lim];
end
line([0 0], ylim,'LineStyle', ':', 'Color', 'k');
line([-0.5 0.5], [0 0],'LineStyle', ':', 'Color', 'k');
xlabel('time (s)');
title(tasklabel);
axis tight
ylabel(signaltoplot);

%% plot individual stimuli: soc

tasklabel = 'soc'; 
signaltoplot = 'broadband';
eltoplot = 'EEG OC_01-REF';

usebootstrap = 'yes'; % yes or no

el = ecog_matchChannels(eltoplot,trials);
temp = squeeze(trials.(signaltoplot).(tasklabel)(el,:,:)); 

% compute STIMULUS averages

nStim = length(category_names);
ga = nan(nStim, length(trials.time)); 
ulim = nan(nStim , length(trials.time));  
llim = nan(nStim , length(trials.time)); 
for stim = 1:nStim 
    categoryInx = find(trials.stimuli.(tasklabel) == stim);
    switch usebootstrap
        case 'yes'
            for tt = 1:length(trials.time)
                [bootstat] = bootstrp(100,@mean,squeeze(temp(tt,categoryInx)));
                ga(stim,tt) = median(bootstat);
                llim(stim,tt) = quantile(sort(bootstat), 0.05);
                ulim(stim,tt) = quantile(sort(bootstat), 0.95);
            end
        otherwise
             ga(stim,:) = mean(temp(:,categoryInx),2);
             llim(stim,:) = ga(stim,:) - std(temp(:,categoryInx),0,2)';
             ulim(stim,:) = ga(stim,:) + std(temp(:,categoryInx),0,2)';
    end 
end

% plot
figure('Name', [tasklabel ' ' signaltoplot ' ' trials.label{el}]); hold on
lim = max(max(abs(ga)));
for stim = 1:nStim
    subplot(ceil(sqrt(nStim)),ceil(sqrt(nStim)),stim); hold on;

    h = ciplot(llim(stim,:),ulim(stim,:),trials.time,'k',0.1);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    plot(trials.time, ga(stim,:), 'k', 'LineWidth',2);
    switch signaltoplot
        case 'broadband'
            ylim = [0 lim+0.2*lim];
        case 'evoked'
            ylim = [-lim-0.2*lim lim+0.2*lim];
    end
    line([0 0], ylim,'LineStyle', ':', 'Color', 'k');
    line([trials.time(1) trials.time(end)], [0 0],'LineStyle', ':', 'Color', 'k');
    set(gca, 'YLim', ylim);
    xlabel('time (s)');
    title(category_names{stim})
end

% plot temporal stimuli 

nStim = length(category_names);
catIndex = 25:nStim;
%colors1 = cool(length(catIndex)/2);
%colors2 = hot(length(catIndex)/2);
%colors = [colors1;colors2];
colors = copper(length(catIndex)/2);

% % in one plot
% figure;hold on
% for stim = 1:length(catIndex)/2
%     h = ciplot(llim(catIndex(stim),:),ulim(catIndex(stim),:),trials.time,colors(stim,:),0.1);
%     h.Annotation.LegendInformation.IconDisplayStyle = 'off';
% end
% for stim = 1:length(catIndex)/2
%     h = ciplot(llim(catIndex(stim)+length(catIndex)/2,:),ulim(catIndex(stim)+length(catIndex)/2,:),trials.time,colors(stim,:),0.1);
%     h.Annotation.LegendInformation.IconDisplayStyle = 'off';
% end
% 
% for stim = 1:length(catIndex)/2
%     plot(trials.time, ga(catIndex(stim),:), 'Color', colors(stim,:), 'LineWidth',2);
% end
% for stim = 1:length(catIndex)/2
%     plot(trials.time, ga(catIndex(stim)+length(catIndex)/2,:), 'Color', colors(stim,:), 'LineStyle', ':', 'LineWidth',2);
% end
% 
% %ylim = [-lim-0.2*lim lim+0.2*lim];
% ylim = [-0.1 lim+0.2*lim];
% line([0 0], ylim,'LineStyle', ':', 'Color', 'k');
%     line([trials.time(1) trials.time(end)], [0 0],'LineStyle', ':', 'Color', 'k');
% xlabel('time (s)');
% axis tight
% legend(category_names(catIndex));
% title([eltoplot ' ' signaltoplot])
% ylabel(signaltoplot);

% in two plots
figure;hold on
for stim = 1:length(catIndex)/2
    h = ciplot(llim(catIndex(stim),:),ulim(catIndex(stim),:),trials.time,colors(stim,:),0.1);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
for stim = 1:length(catIndex)/2
    plot(trials.time, ga(catIndex(stim),:), 'Color', colors(stim,:), 'LineWidth',2);
end

%ylim = [-lim-0.2*lim lim+0.2*lim];
ylim = [-0.1 lim+0.2*lim];
line([0 0], ylim,'LineStyle', ':', 'Color', 'k');
line([trials.time(1) trials.time(end)], [0 0],'LineStyle', ':', 'Color', 'k');
xlabel('time (s)');
axis tight
legend(category_names(catIndex(1:6)));
title([eltoplot ' ' signaltoplot])
ylabel(signaltoplot);


figure; hold on
for stim = 1:length(catIndex)/2
    h = ciplot(llim(catIndex(stim)+length(catIndex)/2,:),ulim(catIndex(stim)+length(catIndex)/2,:),trials.time,colors(stim,:),0.1);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end

for stim = 1:length(catIndex)/2
    plot(trials.time, ga(catIndex(stim)+length(catIndex)/2,:), 'Color', colors(stim,:), 'LineWidth',2);
end

%ylim = [-lim-0.2*lim lim+0.2*lim];
ylim = [-0.1 lim+0.2*lim];
line([0 0], ylim,'LineStyle', ':', 'Color', 'k');
    line([trials.time(1) trials.time(end)], [0 0],'LineStyle', ':', 'Color', 'k');
xlabel('time (s)');
axis tight
legend(category_names(catIndex(7:end)));
title([eltoplot ' ' signaltoplot])
ylabel(signaltoplot);

%% compute SNR per electrode 

SNRepoch = [0.1 0.4]; % e.g. 80 to 400 ms
signaltoplot = 'evoked';
tasklabel = 'soc';

% load NY648_Elec2Node; % elec names from SOM txt file are not in same
% format as in data.label which is annoying; manually labeling electrodes
% with a node match in retinotopic atlas regions for now...
eltocolor1 = {'MO_01', 'MO_02', 'MO_03', 'MO_04', 'G_03', 'G_04', 'G_09', 'G_10', 'G_11', 'G_12', 'G_17', 'G_18', ...
             'DLPA_07', 'DLPA_08', 'DPSL_05'};
eltocolor2 = {'IO_01', 'IO_02', 'IO_03', 'IO_04'};

[~,t1] = min(abs(trials.time-SNRepoch(1)));
[~,t2] = min(abs(trials.time-SNRepoch(2)));
%temp = squeeze(sum(abs(trials.(signaltoplot).(tasklabel)(:,t1:t2,:)),2));
temp = squeeze(mean(trials.(signaltoplot).(tasklabel)(:,t1:t2,:),2));
SNR1 = abs(mean(temp,2));
SNR2 = SNR1./mean(std(trials.(signaltoplot).(tasklabel)(:,1:52,:),0,3),2); % divide by prestim std
SNR3 = SNR1./std(temp,0,2); % divide by std across trials

eltoplot = good_channels;
tempdata = data;
tempdata.label = tempdata.label(eltoplot);

% xaxislabels
for ee = 1:length(eltoplot)
    elnames{ee} = data.label{eltoplot(ee)}(4:end-4);
end

SNR = SNR2;
figure('Name', ['SNR per electrode: ' tasklabel]);hold on

el1 = ecog_matchchannels(eltocolor1,tempdata);
bar(1:length(eltoplot),SNR(eltoplot,:));
bar(el1,SNR(eltoplot(el1)),'FaceColor', 'r');
el2 = ecog_matchchannels(eltocolor2,tempdata);
bar(el2,SNR(eltoplot(el2)),'FaceColor', 'g');

legend({'non-matched electrodes', 'matched electrodes (wang atlas)', 'ventral'});
set(gca,'XTick', 1:126, 'XTickLabel', elnames, 'XTickLabelRotation', 90);
axis tight; ylabel(signaltoplot);
title(['interval ' num2str(SNRepoch(1)*1000) ' to ' num2str(SNRepoch(2)*1000) ' ms: abs(mean)/std prestim']);

%% average and compute an SNR measure per stim category

tasklabel = 'soc'; 
signaltoplot = 'broadband';
eltoplot = 'OC_01';

usebootstrap = 'no'; % yes or no

el = ecog_matchchannels(eltoplot,data);
temp = squeeze(trials.(signaltoplot).(tasklabel)(el,:,:)); 

% compute CATEGORY averages (combining all FACES, LETTERS, SCENES, and ONE_PULSE stimuli)
catList = [1:12 13 17 21 25 31:36 37];
nCat = length(catList);
ga = nan(nCat,length(trials.time)); 
ulim = nan(nCat, length(trials.time));  
llim = nan(nCat, length(trials.time)); 

for cCat = 1:nCat-1
    categoryInx = find(ismember(stimuli.(tasklabel), catList(cCat):catList(cCat+1)-1));
    disp(length(categoryInx));
    switch usebootstrap
        case 'yes'
            for tt = 1:length(trials.time)
                [bootstat] = bootstrp(100,@mean,squeeze(temp(tt,categoryInx)));
                ga(cCat,tt) = median(bootstat);
                llim(cCat,tt) = quantile(sort(bootstat), 0.05);
                ulim(cCat,tt) = quantile(sort(bootstat), 0.95);
            end
        otherwise
             ga(cCat,:) = mean(temp(:,categoryInx),2);
             llim(cCat,:) = ga(cCat,:) - std(temp(:,categoryInx),0,2)';
             ulim(cCat,:) = ga(cCat,:) + std(temp(:,categoryInx),0,2)';
    end 
end

% plot
figure('Name', [tasklabel ' ' signaltoplot ' ' data.label{el}]); hold on
lim = max(max(abs(temp)));
for cCat = 1:nCat-1
    subplot(ceil(sqrt(nCat)),ceil(sqrt(nCat)),cCat); hold on;

    h = ciplot(llim(cCat,:),ulim(cCat,:),trials.time,'k',0.1);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    plot(trials.time, ga(cCat,:), 'k', 'LineWidth',2);
    ylim = [-lim-0.2*lim lim+0.2*lim];
    set(gca, 'YLim', ylim);
    line([0 0], ylim,'LineStyle', ':', 'Color', 'k');
    line([-0.1 0.5], [0 0],'LineStyle', ':', 'Color', 'k');
    xlabel('time (s)');
    axis tight
    title(category_names{catList(cCat)})
end

%% compute onsets per electrode 
% bootstrap and test for zero-crossing?

%% spectra by category
% is there a gamma response in ventral cortex?
% for ventral/lateral electrodes, is there stimulus selectivity (e.g.
% face/scene)?



