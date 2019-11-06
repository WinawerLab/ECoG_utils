
% SCRIPT DESCRIPTION %
% This script Takes BAIR data from NYU School of Medicine, gets onsets,
% writes out separate runs for each tasks, including tsv event files, and
% required BIDS metadata (coordsystem json and electrodes and channels tsv
% files). It is meant to be run cell-by-cell because some manual inputs are
% required for trigger channel selection and identifying noisy channels.

% Check whether we have the ECoG_utils repository on the path
if ~exist('createBIDS_ieeg_json_nyuSOM.m')
    tbUse ECoG_utils;
end

makePlot = 0;

%% Define paths and BIDS specs %%

% Input paths specs
patientID   = 705; % Specify patient's raw folder name here
RawDataDir  = '/Volumes/server/Projects/BAIR/Data/Raw/ECoG/';
BIDSDataDir = '/Volumes/server/Projects/BAIR/Data/BIDS/';

% BIDS specs: assuming defaults for a first session, full visual set:
projectName = 'motor';
sub_label   = ['som' num2str(patientID)]; 
ses_label   = 'nyuecog01';
ses_labelt1 = 'som3t01';
task_label  = {'boldhand',...
               'fingermappingleft', ...
               'gestures', ...
               'boldsat1', ...
               'boldsat2', ...            
               'boldsat3', ... 
               'boldsat4'}; %   
run_label = {'01','01','01','01','01','01','01'}; 

%% DEFINE PATHS AND DATA

% Define paths
[dataReadDir, dataWriteDir, stimWriteDir, T1WriteDir, preprocDir] = bidsconvert_getpaths(patientID, RawDataDir, ...
    BIDSDataDir, projectName, sub_label, ses_label, ses_labelt1);

% Read ECoG data
[data, hdr] = bidsconvert_readecogdata(dataReadDir, ses_label);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% START OF MANUAL SECTION %%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run the cells below by adding several manual inputs

%% DATA TRIM (PATIENT- and SESSION-SPECIFIC)

% Data not trimmed for this patient, because of specific time points used
% to remove certain triggers (see below)

%% CHANNEL IDENTIFICATION

% Manually click through each channel to identify to trigger channel, as
% well as bad channels (specify in next cell)

% Define time axis (in seconds). First time point = 0 (this is assumed by
% the function we used to detect triggers below, and also in fieldtrip).
t = ((0:hdr.nSamples-1)/hdr.Fs); 

% Define the trigger channel name (probably a 'DC' channel, see hdr.label).
triggerChannelName = 'DC1';
triggerChannel = find(strcmp(triggerChannelName,hdr.label));

% Plot the raw voltage time course of each channel
if makePlot
    for cChan = 1:1:size(data,1) 
        figure;plot(t,data(cChan,:)); 
        title([num2str(cChan) ': ' hdr.label{cChan}]);
        xlabel('Time (s)'); ylabel('Raw amplitude (microV)'); set(gca,'fontsize',16); 
        waitforbuttonpress; close; 
    end 
end

%% WRITE DOWN THE FOLLOWING

% This is a list of to be excluded channels that are NOT from the HD grid
exclude_inx = [24 33 50 65 97 104 106 148 173 180 200]; % e.g. 6 12 13 87

% Specify reasons for marked as bad, e.g. spikes, elipeptic,
% outlierspectrum, lowfreqdrift
BADCHANNELS_MANUALTABLE = [num2cell(exclude_inx)' repmat({'spikes'}, [length(exclude_inx) 1])]; 

%% CHECK THE CHANNEL SELECTIONS and save figures

%% [1] Look at spectra without excluding ANY channels:
badChannels = [];

% Generate spectral plot; check command window output for outliers; 
if makePlot 
    inx_notEEGchans = [find(contains(hdr.chantype, 'ecg')); find(contains(hdr.label, {'DC', 'SG', 'Pleth', 'PR', 'OSAT', 'TRIG'}))];
    chansToPlot = setdiff(1:length(hdr.label),[inx_notEEGchans; badChannels]);
    [outliers] = ecog_plotChannelSpectra(data, chansToPlot, hdr); title('All channels spectra');
    saveas(gcf, fullfile(preprocDir, 'figures', 'bidsconversion', sprintf('%s-%s-spectra_allchannels',sub_label, ses_label)), 'epsc');
end

% Plot the timeseries of all channels
if makePlot 
    figure('Name', 'All channels time course');
    plot(t,data(chansToPlot,:)); xlabel('Time (s)'); ylabel('Raw amplitude (microV)'); title('All channels timecourse'); set(gca,'fontsize',16); 
    % Get trigger time points from data file
    [trigger_onsets] = bidsconvert_findtriggers(data, hdr, triggerChannel, 0);
    hold on; plot(trigger_onsets, ones(length(trigger_onsets),1),'k.','MarkerSize', 25, 'LineStyle','none');
    saveas(gcf, fullfile(preprocDir, 'figures', 'bidsconversion', sprintf('%s-%s-timecourse-allchannels',sub_label, ses_label)), 'epsc');
end

% NOTE: Outliers (identified as channels with mean power that is more that
% two standard deviations above or below the average across channels) should
% not be used to automatically identify bad channels, because channels with
% strong activation can have higher power on average! Instead, look at
% the time courses of those channels again to see what makes them stand
% out, and if you think more or less channels need to be excluded, adapt
% the previous cell and run this one again, until satisfied with selection.

if makePlot 
    for cChan = 1:length(outliers)
        figure('Name', sprintf('Outlier %d', cChan)); 
        plot(t, data(outliers(cChan),:)); 
        title([num2str(outliers(cChan)) ': ' hdr.label{outliers(cChan)}]); 
        xlabel('Time (s)'); ylabel('Raw amplitude (microV)'); set(gca,'fontsize',16); 
    end 
end

%% [2] Look at spectra with bad channels excluded:
badChannels = cell2mat(BADCHANNELS_MANUALTABLE(:,1));
badChannelsDescriptions = BADCHANNELS_MANUALTABLE(:,2);

% Generate spectral plot; check command window output for outliers; 
if makePlot 
    inx_notEEGchans = [find(contains(hdr.chantype, 'ecg')); find(contains(hdr.label, {'DC', 'SG', 'Pleth', 'PR', 'OSAT', 'TRIG'}))];
    chansToPlot = setdiff(1:length(hdr.label),[inx_notEEGchans; badChannels]);
    [outliers] = ecog_plotChannelSpectra(data, chansToPlot, hdr); title('Good channels spectra');
    saveas(gcf, fullfile(preprocDir, 'figures', 'bidsconversion', sprintf('%s-%s-spectra_goodchannels',sub_label, ses_label)), 'epsc');
end

% Plot the timeseries of all the good channels, are all noisy moments gone?
if makePlot 
    figure('Name', 'Good channels time course');
    plot(t,data(chansToPlot,:)); xlabel('Time (s)'); ylabel('Raw amplitude (microV)'); title('all good channels'); set(gca,'fontsize',16); 
    % Get trigger time points from data file
    [trigger_onsets] = bidsconvert_findtriggers(data, hdr, triggerChannel, 0);
    hold on; plot(trigger_onsets, ones(length(trigger_onsets),1),'k.','MarkerSize', 25, 'LineStyle','none');
    saveas(gcf, fullfile(preprocDir, 'figures', 'bidsconversion', sprintf('%s-%s-timecourse-goodchannels',sub_label, ses_label)), 'epsc');
end

% NOTE: Outliers (identified as channels with mean power that is more that
% two standard deviations above or below the average across channels) should
% not be used to automatically identify bad channels, because channels with
% strong activation can have higher power on average! Instead, look at
% the time courses of those channels again to see what makes them stand
% out, and if you think more or less channels need to be excluded, adapt
% the previous cell and run this one again, until satisfied with selection.

if makePlot 
    for cChan = 1:length(outliers)
        figure('Name', sprintf('Outlier %d', cChan)); 
        plot(t, data(outliers(cChan),:)); 
        title([num2str(outliers(cChan)) ': ' hdr.label{outliers(cChan)}]); 
        xlabel('Time (s)'); ylabel('Raw amplitude (microV)'); set(gca,'fontsize',16); 
    end 
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% END OF MANUAL SECTION %%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% From here on, everything should run automatically:

% AUTOMATED EXTRACTION %%

%% PATIENT-SPECIFIC HACK ON TRIGGERS

peakOpts.minPeakHeight = 0.3;
peakOpts.minPeakDistance = 0.1;
    
% Get trigger time points from data file
[trigger_onsets] = bidsconvert_findtriggers(data, hdr, triggerChannel, makePlot, [], peakOpts);

if makePlot
    saveas(gcf, fullfile(preprocDir, 'figures', 'bidsconversion', sprintf('%s-%s-triggers_found_all',sub_label, ses_label)), 'epsc');
end

% REMOVE triggers that belong to the practice runs:
% determined by visual inspection that bold hand started after 800 s:
temp_trigger_onsets1 = trigger_onsets(trigger_onsets > 800); 

% REMOVE triggers that belong to the excluded (broken off) fingermapping
% right run:
temp_trigger_onsets2 = trigger_onsets(trigger_onsets < 1380 | trigger_onsets > 1420); 

temp_trigger_onsets = intersect(temp_trigger_onsets1, temp_trigger_onsets2);
trigger_onsets = temp_trigger_onsets;
 
triggers = data(triggerChannel,:);
triggers = triggers / max(triggers);

if makePlot
    figure('Name', 'Found triggers corrected'); hold on;
    plot(t,triggers);
    plot(trigger_onsets, ones(length(trigger_onsets),1)/4,'r.','MarkerSize', 25, 'LineStyle','none');
    legend({'trigger data', 'trigger onsets'}); xlabel('time (s)'); ylabel('amplitude');
    title('found triggers');
    xlabel('Time (s)'); set(gca, 'FontSize', 12);
    saveas(gcf, fullfile(preprocDir, 'figures', 'bidsconversion', sprintf('%s-%s-triggers_found_corrected',sub_label, ses_label)), 'epsc');
end

%% Continue

% Generate electrode files
[electrode_table, channel_table] = bidsconvert_getelectrodefiles(dataReadDir, hdr, triggerChannel, badChannels, badChannelsDescriptions);

% Read in stimulus files
[stimData] = bidsconvert_matchstimulusfiles(dataReadDir, patientID, ses_label, task_label, run_label, trigger_onsets, 1);
if makePlot
    saveas(gcf, fullfile(preprocDir, 'figures', 'bidsconversion', sprintf('%s-%s-triggers_requested',sub_label, ses_label)), 'epsc');
end

% WRITING OF FILES %%%

% Write session files
bidsconvert_writesessionfiles(dataReadDir, dataWriteDir, T1WriteDir, sub_label, ses_label, ses_labelt1, electrode_table)

% Write run files
bidsconvert_writerunfiles(dataWriteDir, stimWriteDir, sub_label, ses_label, task_label, run_label, ...
    data, hdr, stimData, channel_table, trigger_onsets)



