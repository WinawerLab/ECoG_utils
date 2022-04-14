function [channel_table] = createBIDS_ieeg_channels_tsv_nyuSOM(n)

% Adapted from BIDS_starter_kit createBIDS_ieeg_json template
% Intended to be run to get default options for NYU SOM recordings
%
% IGroen, 2018

%% Template Matlab script to create an BIDS compatible channels.tsv file
% This example lists all required and optional fields.
% When adding additional metadata please use CamelCase 
%
% DHermes, 2017
% modified RG 201809

%% required columns
name = repmat('n/a', n, 1); % Label of the channel, only contains letters and numbers. 
% The label must correspond to _electrodes.tsv name and all ieeg type
% channels are required to have a position. The reference channel name MUST
% be provided in the reference column.

type = repmat('OTHER', n, 1); % Type of channel, see section 3.3.2 of iEEG BIDS spec.

units = repmat({'microV'}, n, 1); % Physical unit of the value represented in this channel,
% e.g. V for Volt, specified according to the SI unit symbol and possibly 
% prefix (e.g. milliV, microV), see BIDS spec for Units and Prefixes.

low_cutoff = repmat('500', n, 1); %Frequencies used for the low pass filter applied to the 
% channel in Hz. If no low pass filter was applied, use n/a. Note that 
% anti-alias is a low pass filter, specify its frequencies here if applicable.

high_cutoff = repmat('0.15', n, 1); % Frequencies used for the high pass filter applied to 
% the channel in Hz. If no high pass filter applied, use n/a.

%% recommended and optional columns:
sampling_frequency = repmat('unknown', n, 1); % OPTIONAL. Sampling rate of 
% the channel in Hz

notch = repmat('n/a', n, 1); % OPTIONAL. Frequencies used for the notch 
% filter applied to the channel, in Hz. If no notch filter applied, use n/a 

reference = repmat('n/a', n, 1); % Specification of the reference (options: ?bipolar?, 
% ?mastoid?, ?intracranial?, ?ElectrodeName01?) .

group = repmat({'n/a'}, n, 1); % RECOMMENDED Which group of channels 
% (grid/strip/seeg/depth) this channel belongs to. This is relevant because
% one group has one cable-bundle and noise can be shared. This can be a
% name or number. Note that any groups specified in `_electrodes.tsv` must
% match those present here

description = repmat({'n/a'}, n, 1); % OPTIONAL. Brief free-text description 
% of the channel, or other information of interest (e.g. position (e.g.,
% ?left lateral temporal surface?, etc.)

status = repmat({'good'}, n, 1); % OPTIONAL. Data quality observed on the 
% channel (good/bad). A channel is considered bad if its data quality is
% compromised by excessive noise. Description of noise type SHOULD be
% provided in [status_description].

status_description = repmat({'n/a'}, n, 1); % OPTIONAL. Freeform text description of noise 
% or artifact affecting data quality on the channel. It is meant to explain
% why the channel was declared bad in [status].

%% write
channel_table = table(name,type,units,low_cutoff,high_cutoff, sampling_frequency, notch, reference, ...
    group,description,status,status_description);
