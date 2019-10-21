
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

%% Define paths and BIDS specs %%

% Input paths specs
patientID   = 748; % Specify patient's raw folder name here
RawDataDir  = '/Volumes/server/Projects/BAIR/Data/Raw/ECoG/';
BIDSDataDir = '/Volumes/server/Projects/BAIR/Data/BIDS/';

% BIDS specs: assuming defaults for a first session, full visual set:
projectName = 'visual';
sub_label   = ['som' num2str(patientID)]; 
ses_label   = 'nyuecog02';
ses_labelt1 = 'som3t01';
task_label  = {'prf',...
               'prf', ...
               'hrfpattern', ...
               'temporalpattern', ...            
               'temporalpattern', ... 
               'spatialpattern', ...
               'spatialpattern', ...    
               'spatialobject', ...
               'spatialobject', ...
               'temporalpattern', ...            
               'temporalpattern', ...       
               'spatialpattern', ...
               'spatialpattern', ...
               'spatialobject', ...
               'spatialobject'
              };              
run_label = {'01','02','01','01','02','01','02','01','02','03','04','03','04','03','04'};
% NOTE: task and run labels should be noted in the order they were run!

% Make plots?
makePlot = 1;
% NOTE: Figures will be saved into
% derivatives/preprocessed/sub-label/ses-label/figures/bidsconversion

%% DEFINE PATHS AND DATA

% Define paths
[dataReadDir, dataWriteDir, stimWriteDir, T1WriteDir, preprocDir] = bidsconvert_getpaths(patientID, RawDataDir, ...
    BIDSDataDir, projectName, sub_label, ses_label, ses_labelt1);

% Read ECoG data
[rawdata, hdr] = bidsconvert_readecogdata(dataReadDir, ses_label);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% START OF MANUAL SECTION %%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run the cells below by adding several manual inputs

%% DATA TRIM (PATIENT- and SESSION-SPECIFIC)

% Define the trigger channel name (probably a 'DC' channel, see hdr.label).
triggerChannelName = 'DC1';
triggerChannel = find(strcmp(triggerChannelName,hdr.label));
figure;plot(rawdata(triggerChannel,:)); 
title([num2str(triggerChannel) ': ' hdr.label{triggerChannel}]);
        
run_start = [140000]; % Manually determined from plot of triggerchannel 
run_end   = [992000]; 

% Clip the data
data = rawdata(:,run_start:run_end);
hdr.nSamples = size(data,2);

% Check if we have all the triggers we want
figure;plot(data(triggerChannel,:)); 
title([num2str(triggerChannel) ': ' hdr.label{triggerChannel}]);

% PATIENTSPECIFIC HACK %%

% For 748, there's a mismatch for one set of electrodes that are labeled
% 'DMP' in the data, but 'DMPT' in the electrode coordinates provided by
% SoM. Hack by overwriting names in the hdr:
INX = find(contains(hdr.label, 'DPM')); 
for ii = 1:length(INX)
    oldlabel = hdr.label{INX(ii)};
    newlabel = [oldlabel(1:3) 'T' oldlabel(4:end)];
    hdr.label{INX(ii)} = newlabel;
end

%% CHANNEL IDENTIFICATION

% Manually click through each channel to identify to trigger channel, as
% well as bad channels (specify in next cell)

% Define time axis (in seconds). First time point = 0 (this is assumed by
% the function we used to detect triggers below, and also in fieldtrip).
t = ((0:hdr.nSamples-1)/hdr.Fs); 

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

% This is a list of to be excluded channels from CAR; will be labeled as
% 'bad' in the channels tsv file
exclude_inx = [96 126 144 189 190 191 207 129 248 251]; 

% Specify reasons for marked as bad, e.g. spikes, elipeptic,
% outlierspectrum, lowfreqdrift
%BADCHANNELS_MANUALTABLE = [num2cell(exclude_inx)' repmat({'spikes'}, [length(exclude_inx) 1])]; 
BADCHANNELS_MANUALTABLE = [num2cell(exclude_inx)' {'saturated'; 'saturated'; 'spikes'; 'spikes'; 'spikes'; 'spikes'; 'spikes'; 'spikes'; 'spikes'; 'spikes'}];

%% CHECK THE CHANNEL SELECTIONS and save figures

%% [1] Look at spectra without excluding ANY channels:
badChannels = [];

% Generate spectral plot; check command window output for outliers; 
if makePlot 
    inx_notEEGchans = find(contains(hdr.label, {'DC', 'SG', 'Pleth', 'PR', 'OSAT', 'TRIG', 'EKG', 'ECG'}));
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
    inx_notEEGchans = find(contains(hdr.label, {'DC', 'SG', 'Pleth', 'PR', 'OSAT', 'TRIG', 'EKG', 'ECG'}));
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

% Get trigger time points from data file
[trigger_onsets] = bidsconvert_findtriggers(data, hdr, triggerChannel, makePlot);
if makePlot
    saveas(gcf, fullfile(preprocDir, 'figures', 'bidsconversion', sprintf('%s-%s-triggers_found',sub_label, ses_label)), 'epsc');
end

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



