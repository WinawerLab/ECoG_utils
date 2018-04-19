% tbUse('ECoG_utils');

%%%% DEFINE PATHS TO DATA, ANALYSIS DIRECTORIES %%%%

anaDir = '/Volumes/server/Projects/BAIR/Analyses/';
projectName = 'visual';

sub_label = 'beilen'; 
ses_label = 'day03';

%% Load preprocessed data

dataToLoad = fullfile(anaDir, projectName, ['sub-' sub_label], ...
    ['sub-' sub_label '_' 'ses-' ses_label '_epoched']);

load(dataToLoad);

%% Plot trial averages

% According to sub-beilen_acq-clinicalprojectedregions_electrodes.tsv: 
% Electrodes with lateral-occipital coverage: 
eltomatch = {'BT12'}; % 'OT07'
% Electrodes with inferior-parietal coverage: 
%eltomatch = {'BT16', 'OT08', 'OT16'};

signaltoplot = 'broadband';
el = find(strcmp(eltomatch{1}, trials.label));

figure; hold on
title(eltomatch);
plot(trials.time, median(trials.(signaltoplot).bairhrfpattern(el,:,:),3),'b', 'LineWidth', 2)
plot(trials.time, median(trials.(signaltoplot).bairprf(el,:,:),3),'r','LineWidth', 2)
plot(trials.time, median(trials.(signaltoplot).bairspatialobject(el,:,:),3),'m', 'LineWidth', 2)
plot(trials.time, median(trials.(signaltoplot).bairspatialpattern(el,:,:),3),'c', 'LineWidth', 2)
plot(trials.time, median(trials.(signaltoplot).bairtemporalpattern(el,:,:),3),'g', 'LineWidth', 2)
legend(trials.tasks);

% add stim onset and zero lines
ylim = get(gca, 'ylim');
l1 = line([0 0], [ylim],'LineStyle', ':', 'Color', 'k'); l1.Annotation.LegendInformation.IconDisplayStyle = 'off';
l2 = line([-0.5 1.5], [0 0],'LineStyle', ':', 'Color', 'k'); l2.Annotation.LegendInformation.IconDisplayStyle = 'off';

xlabel('time (s)');
switch signaltoplot
    case 'evoked'
        ylabel('evoked amplitude')
    case 'broadband'
        ylabel('broadband (60-120 Hz)');
end

%% plot single trials: hrf

eltomatch = {'OT15'}; % 'OT07'
signaltoplot = 'broadband';
el = find(strcmp(eltomatch{1}, trials.label));

nTrials = size(trials.(signaltoplot).bairhrfpattern,3);
figure;hold on
for ii = 1:nTrials
    subplot(ceil(sqrt(nTrials)), ceil(sqrt(nTrials)), ii); hold on;
    plot(trials.time, trials.(signaltoplot).bairhrfpattern(el,:,ii),'k', 'LineWidth', 1)
    % add stim onset and zero lines
    ylim = [-2 40];
    set(gca, 'ylim', ylim);
    l1 = line([0 0], ylim,'LineStyle', ':', 'Color', 'k'); l1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    l2 = line([-0.5 1.5], [0 0],'LineStyle', ':', 'Color', 'k'); l2.Annotation.LegendInformation.IconDisplayStyle = 'off';
end

figure;hold on
plot(trials.time, squeeze(trials.(signaltoplot).bairhrfpattern(el,:,:)),'k', 'LineWidth', 1)
ylim = [-2 40];
set(gca, 'ylim', ylim);
l1 = line([0 0], ylim,'LineStyle', ':', 'Color', 'k'); l1.Annotation.LegendInformation.IconDisplayStyle = 'off';
l2 = line([-0.5 1.5], [0 0],'LineStyle', ':', 'Color', 'k'); l2.Annotation.LegendInformation.IconDisplayStyle = 'off';

%% plot category averages: 

