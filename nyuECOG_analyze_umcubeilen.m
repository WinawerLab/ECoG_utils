tbUse('ECoG_utils');

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
eltomatch = {'OT15'}; %, 'OT07'
% Electrodes with inferior-parietal coverage: 
%eltomatch = {'BT16', 'OT08', 'OT16'};

signaltoplot = 'evoked';
el = find(strcmp(eltomatch{1}, trials.label));

figure; hold on
plot(trials.time, mean(trials.(signaltoplot).bairhrfpattern(el,:,:),3),'b', 'LineWidth', 2)
plot(trials.time, mean(trials.(signaltoplot).bairprf(el,:,:),3),'r','LineWidth', 2)
plot(trials.time, mean(trials.(signaltoplot).bairspatialobject(el,:,:),3),'m', 'LineWidth', 2)
plot(trials.time, mean(trials.(signaltoplot).bairspatialpattern(el,:,:),3),'c', 'LineWidth', 2)
plot(trials.time, mean(trials.(signaltoplot).bairtemporalpattern(el,:,:),3),'g', 'LineWidth', 2)
legend(allTasks);

% add stim onset and zero lines
l = line([0 0], ylim,'LineStyle', ':', 'Color', 'k'); l.Annotation.LegendInformation.IconDisplayStyle = 'off';
l = line([-0.1 0.5], [0 0],'LineStyle', ':', 'Color', 'k'); l.Annotation.LegendInformation.IconDisplayStyle = 'off';

xlabel('time (s)');
ylabel('broadband (60-120 Hz)');

% TO DO 
% - adapt electrode_to_node to read umcu data, extract electrode locations
% - figure out categories for spatiotemporal runs (refer back to original .tsv files?
% - plot temporal stimuli with ciplot


