
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
patientID   = 723;
RawDataDir  = '/Volumes/server/Projects/BAIR/Data/Raw/ECoG/';
BIDSDataDir = '/Volumes/server/Projects/BAIR/Data/BIDS/';

% BIDS specs
projectName = 'visual';
sub_label   = ['som' num2str(patientID)]; 
ses_label   = 'nyuecog02';
ses_labelt1 = 'som3t01';
task_label  = {'prf',...
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
               'hrfpattern', ...
              };              
run_label = {'01','02','01','02','01','02','01','02','03','04','03','04','03','04', '02','01'};

% Make plots?
makePlots = 1;

%% Define paths and read data

[dataReadDir, dataWriteDir, stimWriteDir, T1WriteDir] = bidsconvert_getpaths(patientID, RawDataDir, ...
    BIDSDataDir, projectName, sub_label, ses_label, ses_labelt1);

[data, hdr] = bidsconvert_readecogdata(dataReadDir, ses_label);

% PATIENTSPECIFIC HACK %%
% For 723, there's a mismatch for one set of electrodes that are labeled
% 'RDS' in the data, but 'RDSI' in the electrode coordinates provided by
% SoM. Hack by overwriting names in the hdr:
INX = find(contains(hdr.label, 'RDS')); 
for ii = 1:length(INX)
    oldlabel = hdr.label{INX(ii)};
    newlabel = [oldlabel(1:3) 'I' oldlabel(4:end)];
    hdr.label{INX(ii)} = newlabel;
end

%% START OF MANUAL SECTION %%

% Use the snippets of code from this cell to do trigger channel selection,
% and to identify good and bad channels (specify in next cell)

% Define time axis (in seconds). First time point = 0 (this is assumed by
% the function we used to detect triggers below, and also in fieldtrip).
t = ((0:hdr.nSamples-1)/hdr.Fs); 

% Plot the raw voltage time course of each channel
if makePlots, for cChan = size(data,1):-1:1; figure;plot(t,data(cChan,:)); title([num2str(cChan) ': ' hdr.label{cChan}]); waitforbuttonpress; close; end; end

% WRITE DOWN THE FOLLOWING

% Trigger channel: index for the trigger channel (probably one labeled 'DC', see hdr.label)
triggerChannel = find(strcmp('DC1',hdr.label));

% Bad channels (e.g. those with big spikes):
badChannels = find(contains(hdr.label,{'SG', 'Pleth', 'PR', 'OSAT', 'TRIG'}))';

% Indicate the reason why channels were marked as bad
badChannelsDescriptions = repmat({'bad'}, [length(badChannels) 1]);%{'subgenual','subgenual','subgenual','subgenual'};

% All sEEG/ECOG channels not labeled as bad will be labeled good
goodChannels = setdiff(chansToPlot,badChannels)';

% CHECK THE CHANNEL SELECTIONS

% Generate spectral plot; check command window output for outliers; 
inx_notEEGchans = find(contains(hdr.chantype, 'ecg'));
inx_DCchans = find(contains(hdr.label, 'DC'));
chansToPlot = setdiff(1:length(hdr.label),[inx_notEEGchans; inx_DCchans; badChannels']);
[outliers] = ecog_plotChannelSpectra(data, chansToPlot,hdr);

% NOTE: outliers (identified as channels with mean power that is more that
% two standard deviations above or below the average across channels) should
% not be used to automatically identify bad channels, because channels with
% strong activation can have higher power on average! Instead, look at
% the time courses of those channels again to see what makes them stand out:

if makePlots, for cChan = 1:length(outliers);figure; plot(t, data(outliers(cChan),:)); title([num2str(outliers(cChan)) ': ' hdr.label{outliers(cChan)}]); end; end

% Check the powerplot of all the good channels, no leftover outliers?
if makePlots, ecog_plotChannelSpectra(data, goodChannels, hdr); end

% Check the timeseries of all the good channels, no noisy moments?
if makePlots, figure;plot(t,data(goodChannels,:)); xlabel('Time (s)'); ylabel('Amplitude'); set(gca,'fontsize',16); end

%% END OF MANUAL SECTION %%

% From here on, everything should run automatically:

%% AUTOMATED EXTRACTION %%

% Generate electrode files
[electrode_table, channel_table] = bidsconvert_getelectrodefiles(dataReadDir, hdr, triggerChannel, goodChannels, badChannels, badChannelsDescriptions);

% Get trigger time points from data file
[trigger_onsets] = bidsconvert_findtriggers(data, hdr, triggerChannel, makePlots);

% Read in stimulus files
[stimData] = bidsconvert_readstimulusfiles(dataReadDir, patientID, ses_label, run_label, trigger_onsets, makePlots);

%% WRITING %%%

% Write session files
bidsconvert_writesessionfiles(dataReadDir, dataWriteDir, T1WriteDir, sub_label, ses_label, ses_labelt1, electrode_table)

% Write run files
bidsconvert_writerunfiles(dataWriteDir, stimWriteDir, sub_label, ses_label, task_label, run_label, ...
    data, hdr, stimData, channel_table, trigger_onsets)



