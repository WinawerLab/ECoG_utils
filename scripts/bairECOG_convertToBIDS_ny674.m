
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
patientID   = 674; % Specify patient's raw folder name here
RawDataDir  = '/Volumes/server/Projects/BAIR/Data/Raw/ECoG/';
BIDSDataDir = '/Volumes/server/Projects/BAIR/Data/BIDS/';

% BIDS specs: assuming defaults for a first session, full visual set:
projectName = 'visual';
sub_label   = ['som' num2str(patientID)]; 
ses_label   = 'nyuecog01';
ses_labelt1 = 'som3t01';
acq_label   = 'clinical';
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
               'prf',...
               'prf'
              };              
run_label = {'01','01','02','01','02','01','02','01','02','03','04','03','04','03','04', '03', '04'};

%% DEFINE PATHS AND DATA

% Define paths
[dataReadDir, dataWriteDir, stimWriteDir, T1WriteDir, preprocDir] = bidsconvert_getpaths(patientID, RawDataDir, ...
    BIDSDataDir, projectName, sub_label, ses_label, ses_labelt1);

% PATIENT SPECIFIC NOTE:
% The data for this patient was split up across three different .edf files,
% with two files containing data from the same set of electrodes but split
% it two different time chunks, and the third file containing another set
% of electrodes with the entire time course. The two different sets of
% electrodes data were also segmented out from the continuous recording
% starting at different time points so we need to use custom code to read
% in the data and align the two electrode sets in time.

[onsets, onsets_run, data, hdr] = NY674_getOnsetsHeaderData();

% Make sure there is some padding in between runs
Fs = hdr{1}.samples(1);
nSamplesPre = 3*Fs;
nSamplesPost = 3*Fs;

% To concatenate data in electrode dimension, we need to start and stop at
% same triggers:

data_new = cell(size(data));
trigger_onsets = [];
totalRunLengthInSamples = 0;
for ii = 1:length(onsets_run{1})
    % determine how long the run should last
    runLengthInSamples = [];
    for dd = 1:2
        firsttrig = onsets_run{dd}{ii}(1)-nSamplesPre;
        lasttrig = onsets_run{dd}{ii}(end)+nSamplesPost;
        runLengthInSamples(dd) = lasttrig-firsttrig;
    end
    runLengthInSamples = max(runLengthInSamples);
    % segment each run according to the same runlength
    for dd = 1:2
        firsttrig = onsets_run{dd}{ii}(1)-nSamplesPre;
        data_new{dd} = [data_new{dd} data{dd}(:,firsttrig:firsttrig+runLengthInSamples-1)];
    end
    % recompute the onsets 
    onsets_new = (onsets_run{2}{ii}-firsttrig)+totalRunLengthInSamples;
    trigger_onsets = [trigger_onsets onsets_new]; % in samples
    totalRunLengthInSamples = totalRunLengthInSamples+runLengthInSamples;
end

% Concatenate
nonOverlappingChannels = ~ismember(hdr{1}.label, hdr{2}.label);
data = [data_new{1}(nonOverlappingChannels,:);data_new{2}];

% Create a fieldtrip compatible header:
oldhdr = hdr; hdr = [];
labels = [oldhdr{1}.label(nonOverlappingChannels) oldhdr{2}.label];

hdr.Fs = Fs;
hdr.nChans = length(labels);
hdr.label = labels';
hdr.nSamples = size(data,2);
hdr.nSamplesPre = 0;
hdr.nTrials = 1;
hdr.orig = [];
hdr.chantype = [];
hdr.chanunit = [oldhdr{1}.units(nonOverlappingChannels) oldhdr{2}.units]';

% Convert trigger_onsets to seconds:
trigger_onsets = trigger_onsets * (1/hdr.Fs);

%% START OF MANUAL SECTION %%

% Define time axis (in seconds). First time point = 0 (this is assumed by
% the function we used to detect triggers below, and also in fieldtrip).
t = ((0:hdr.nSamples-1)/hdr.Fs); 

% Plot the raw voltage time course of each channel
if makePlot
    for cChan = 127:1:size(data,1) 
        figure;plot(t,data(cChan,:)); 
        title([num2str(cChan) ': ' hdr.label{cChan}]);
        xlabel('Time (s)'); ylabel('Raw amplitude (microV)'); set(gca,'fontsize',16); 
        waitforbuttonpress; close; 
    end 
end

%% WRITE DOWN THE FOLLOWING

% Trigger channel name (probably a 'DC' channel, see hdr.label)
triggerChannelName = 'DC10';

% Bad channel numbers (e.g. those with big spikes):
BADCHANNELS_MANUALTABLE = {...
    147, 'spikes';
    150, 'excessivenoise';
    157, 'excessivenoise';
    171, 'outlierspectrum';
    185, 'subgenual';
    186, 'subgenual';
    187, 'subgenual';
    188, 'subgenual';
    189, 'subgenual';
    190, 'subgenual';
    };

% NOTE: there is a lot of shared noise on grid A and B, making it
% impossible to see the raw time courses. Preliminary analysis showed that
% this was resolved through the common average referencing. For now,
% labeling all grid channels as good; even though there is a group of
% electrodes with deviating spectra in (one of?) the grids (see below).
% Check CAR in preprocessing and perhaps try CAR per amplifier (see below).

% NOTE from Adeen (per email, May 15 2018), on when data looks identical
% across channels: "The data is fine but sometimes the EEG techs mess up
% the reference a bit. If you run a common average reference it will all
% clear up. In case that doesn?t work (it did for me) you can do a CAR in
% groups of 64 which would remove shared noise per amplifier."

%% CHECK THE CHANNEL SELECTIONS

% Compare selection above with excluding NO channels (comment/uncomment):
%BADCHANNELS_MANUALTABLE  = {[],[]};

triggerChannel = find(contains(hdr.label,triggerChannelName));
badChannels = cell2mat(BADCHANNELS_MANUALTABLE(:,1));
badChannelsDescriptions = BADCHANNELS_MANUALTABLE(:,2);

% Generate spectral plot; check command window output for outliers; 
if makePlot 
    inx_notEEGchans = find(~contains(hdr.label, {'EEG'}))';
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

% USING TRIGGER ONSETS COMPUTED ABOVE
%[trigger_onsets] = bidsconvert_findtriggers(data, hdr, triggerChannel, makePlot);

if makePlot
    saveas(gcf, fullfile(preprocDir, 'figures', 'bidsconversion', sprintf('%s-%s-triggers_found',sub_label, ses_label)), 'epsc');
end

% Generate electrode files
[electrode_table, channel_table] = bidsconvert_getelectrodefiles(dataReadDir, hdr, triggerChannel, badChannels, badChannelsDescriptions);

% Read in stimulus files
[stimData, triggersAreMatched, runTimes] = bidsconvert_matchstimulusfiles(dataReadDir, patientID, ses_label, task_label, run_label, trigger_onsets, 1);
if makePlot
    saveas(gcf, fullfile(preprocDir, 'figures', 'bidsconversion', sprintf('%s-%s-triggers_requested',sub_label, ses_label)), 'epsc');
end

% WRITING OF FILES %%%

% Write run files
[dataFileNames] = bidsconvert_writerunfiles(dataWriteDir, stimWriteDir, ...
    sub_label, ses_label, task_label, acq_label, run_label, ...
    data, hdr, stimData, channel_table, trigger_onsets);

% Write session files
bidsconvert_writesessionfiles(dataReadDir, dataWriteDir, T1WriteDir, ...
    sub_label, ses_label, acq_label, ses_labelt1, electrode_table, dataFileNames, runTimes);



