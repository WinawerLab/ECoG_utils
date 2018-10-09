tbUse('ECoG_utils');

% SCRIPT DESCRIPTION %

%% [0] Define paths and dataset specs %%

% Dataset specs
projectName = 'visual';
sub_label   = 'som648'; 
ses_label   = 'nyuecog01';

% Input paths specs
dataPth     = '/Volumes/server/Projects/BAIR/Data/BIDS/';

% Output paths specs
dataDir = fullfile(dataPth, projectName, sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label), 'ieeg');
saveDir = fullfile(dataPth, projectName, 'derivatives', 'preprocessed', sprintf('sub-%s', sub_label));


%% Load preprocessed data

dataToLoad = fullfile(anaDir, projectName, ['sub-' sub_label], ...
    ['sub-' sub_label '_' 'ses-' ses_label '_' 'task-' task_label '_epoched']);

load(dataToLoad);

%% Plot trial averages for channel(s) of interest

usebootstrap = 'no'; % yes/no
signaltoplot = 'broadband'; % evoked/broadband

% % % electrodes with coverage in V1-V3:
eltomatch = {'MO_01', 'MO_02', 'MO_03', 'MO_04'};
% % % electrodes with coverage in IPS (surface):
% eltomatch = {'G_03', 'G_04', 'G_09', 'G_10', 'G_11', 'G_12', 'G_17', 'G_18'};
% % % electrodes with coverage in IPS (depth):
% eltomatch = {'DLPA_07', 'DLPA_08', 'DPSL_05'};

% % % electrodes near ventrolateral regions:
% eltomatch = {'IO_01', 'IO_02', 'IO_03', 'IO_04'};
% % % electrodes with good SNR in broadband, but no match with atlas:
% eltomatch = {'G_05', 'G_06', 'G_07', 'G_08', 'G_16', 'G_19'};

% find matching electrode numbers
el = ecog_matchchannels(eltomatch, trials);

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
    line([-0.1 0.5], [0 0],'LineStyle', ':', 'Color', 'k');
    
    % add legends, labels
    legend(data.label(el));
    xlabel('time (s)');
    title(tasklabel);
    axis tight
    ylabel(signaltoplot);
end

%% plot individual trials: hrf

tasklabel = 'hrf'; 
signaltoplot = 'broadband'; 
eltoplot = 'MO_01';

el = ecog_matchchannels(eltoplot,trials);
temp = squeeze(trials.(signaltoplot).(tasklabel)(el,:,:)); 

% plot trials in separate panels
figure('Name', [tasklabel ' ' signaltoplot ' ' trials.label{el}]); hold on
lim = max(max(abs(temp)));
ylim = [-lim-0.1*lim lim+0.1*lim];
ntrial = size(temp,2);

for trial = 1:ntrial
    subplot(ceil(sqrt(ntrial)),ceil(sqrt(ntrial)),trial); hold on;
    
    plot(trials.time, temp(:,trial), 'k', 'LineWidth',2);
    set(gca, 'YLim', ylim);
    line([0 0], ylim,'LineStyle', ':', 'Color', 'k');
    line([-0.1 0.5], [0 0],'LineStyle', ':', 'Color', 'k');
    xlabel('time (s)');
    %ylabel(ylabels{signaltoplot})
    axis tight
    if trial == 1
        title([{'trial number:'} {num2str(trial)}])
    else
        title(num2str(trial));
    end
end

% plot trials on top of one another
figure('Name', [tasklabel ' ' signaltoplot ' ' trials.label{el}]); hold on
plot(trials.time,temp);
set(gca, 'YLim', ylim);
line([0 0], ylim,'LineStyle', ':', 'Color', 'k');
line([-0.5 0.5], [0 0],'LineStyle', ':', 'Color', 'k');
xlabel('time (s)');
title(tasklabel);
axis tight
ylabel(signaltoplot);

%% plot individual stimuli: soc

tasklabel = 'soc'; 
signaltoplot = 'broadband';
eltoplot = 'MO_01';

usebootstrap = 'no'; % yes or no

el = ecog_matchchannels(eltoplot,trials);
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
    ylim = [-lim-0.2*lim lim+0.2*lim];
    set(gca, 'YLim', ylim);
    line([0 0], ylim,'LineStyle', ':', 'Color', 'k');
    line([trials.time(1) trials.time(end)], [0 0],'LineStyle', ':', 'Color', 'k');
    xlabel('time (s)');
    axis tight
    title(category_names{stim})
end

%% plot temporal stimuli 
nStim = length(category_names);

catIndex = 25:nStim;
%colors1 = cool(length(catIndex)/2);
%colors2 = hot(length(catIndex)/2);
%colors = [colors1;colors2];
colors = hsv(length(catIndex)/2);

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
eltoplot = 'G_08';

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