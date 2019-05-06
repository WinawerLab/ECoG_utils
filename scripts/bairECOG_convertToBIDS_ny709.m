
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
patientID   = 709; % Specify patient's raw folder name here
RawDataDir  = '/Volumes/server/Projects/BAIR/Data/Raw/ECoG/';
BIDSDataDir = '/Volumes/server/Projects/BAIR/Data/BIDS/';

% BIDS specs: assuming defaults for a first session, full visual set:
projectName = 'visual';
sub_label   = ['som' num2str(patientID)]; 
ses_label   = 'nyuecog01';
ses_labelt1 = 'som3t01';
task_label  = {'hrfpattern', ...
               'prf',...
               'prf', ...
               'spatialpattern', ...
               'spatialpattern', ...
               'spatialobject', ...
               'spatialobject', ...
               'temporalpattern', ...            
               'temporalpattern', ... 
               'spatialpattern', ...
               'spatialpattern', ...
               'spatialobject', ...
               'spatialobject', ...
               'temporalpattern', ...            
               'temporalpattern', ...
              };              
run_label = {'01','01','02','01','02','01','02','01','02','03','04','03','04','03','04'};

% Make plots?
makePlot = 1;

%% DEFINE PATHS AND DATA

% Define paths
[dataReadDir, dataWriteDir, stimWriteDir, T1WriteDir, preprocDir] = bidsconvert_getpaths(patientID, RawDataDir, ...
    BIDSDataDir, projectName, sub_label, ses_label, ses_labelt1);

% Read ECoG data
[data, hdr] = bidsconvert_readecogdata(dataReadDir, ses_label);

%% START OF MANUAL SECTION %%

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

% Trigger channel name (probably a 'DC' channel, see hdr.label)
triggerChannelName = 'DC1';

% Bad channel numbers (e.g. those with big spikes):
    % First column: bad channel index (number)
        % (all sEEG/ECOG channels not labeled as bad will be labeled good)
    % Second column: indicate reason why marked as bad (optional)
        % (status channels and SG, DC, ECG will be excluded automatically)

% EXAMPLE
BADCHANNELS_MANUALTABLE = {...
    8, 'spikes+lowfreq';
    77, 'spikes';
    88, 'spikes';
    83, 'spikes';
    93, 'spikes';
    94, 'spikes';
    108, 'spikes';
    109, 'spikes';
    150, 'spikes';
    129, 'excessivenoise';
    130, 'excessivenoise';
    131, 'excessivenoise';
    132, 'excessivenoise';
    133, 'excessivenoise';
    134, 'excessivenoise';
    135, 'excessivenoise';
    136, 'excessivenoise';
    137, 'excessivenoise';
    138, 'excessivenoise';
    139, 'excessivenoise';
    140, 'excessivenoise';
    141, 'excessivenoise';
    142, 'excessivenoise';
    143, 'excessivenoise';
    144, 'excessivenoise';
    145, 'excessivenoise';
    146, 'excessivenoise';
    147, 'excessivenoise';
    148, 'excessivenoise';
    149, 'excessivenoise';
    151, 'excessivenoise';
    152, 'excessivenoise';
    153, 'excessivenoise';
    154, 'excessivenoise';
    155, 'excessivenoise';
    156, 'excessivenoise';
    157, 'excessivenoise';
    158, 'excessivenoise';
    159, 'excessivenoise';
    160, 'excessivenoise';
    161, 'excessivenoise';
    162, 'excessivenoise';
    163, 'excessivenoise';
    164, 'excessivenoise';
    165, 'excessivenoise';
    166, 'excessivenoise';
    167, 'excessivenoise';
    168, 'excessivenoise';
    169, 'excessivenoise';
    170, 'excessivenoise';
    171, 'excessivenoise';
    172, 'excessivenoise';
    173, 'excessivenoise';
    174, 'excessivenoise';   
    };

%% CHECK THE CHANNEL SELECTIONS

triggerChannel = find(strcmp(triggerChannelName,hdr.label));
badChannels = cell2mat(BADCHANNELS_MANUALTABLE(:,1));
badChannelsDescriptions = BADCHANNELS_MANUALTABLE(:,2);

% Generate spectral plot; check command window output for outliers; 
if makePlot 
    inx_notEEGchans = [find(contains(hdr.chantype, 'ecg')); find(contains(hdr.label, {'DC', 'SG', 'Pleth', 'PR', 'OSAT', 'TRIG'}))];
    chansToPlot = setdiff(1:length(hdr.label),[inx_notEEGchans; badChannels]);
    [outliers] = ecog_plotChannelSpectra(data, chansToPlot, hdr); title('channel spectra');
    saveas(gcf, fullfile(preprocDir, 'figures', 'bidsconversion', sprintf('%s-%s-goodchannels_spectra',sub_label, ses_label)), 'epsc');
end

% Check the timeseries of all the good channels, no noisy moments?
if makePlot 
    figure('Name', 'Good channels time course');
    plot(t,data(chansToPlot,:)); xlabel('Time (s)'); ylabel('Raw amplitude (microV)'); title('all good channels'); set(gca,'fontsize',16); 
    % Get trigger time points from data file
    [trigger_onsets] = bidsconvert_findtriggers(data, hdr, triggerChannel, 0);
    hold on; plot(trigger_onsets, ones(length(trigger_onsets),1),'k.','MarkerSize', 25, 'LineStyle','none');
    saveas(gcf, fullfile(preprocDir, 'figures', 'bidsconversion', sprintf('%s-%s-goodchannels_timecourse',sub_label, ses_label)), 'epsc');
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

%% END OF MANUAL SECTION %

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



