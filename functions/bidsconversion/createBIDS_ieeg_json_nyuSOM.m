function [ieeg_json, json_options] = createBIDS_ieeg_json_nyuSOM()

% Adapted from BIDS_starter_kit createBIDS_ieeg_json template
% Intended to be run to get default options for NYU SOM recordings
%
% IGroen, 2018

%% Template Matlab script to create an BIDS compatible _ieeg.json file
% This example lists all required and optional fields.
% When adding additional metadata please use CamelCase
%
% DHermes, 2017
% modified RG 201809


%% General fields, shared with MRI BIDS and MEG BIDS:



%% Required fields:

ieeg_json.TaskName = 'n/a'; % Name of the task (for resting state use the rest
% prefix). No two tasks should have the same name. Task label is derived 
% from this field by removing all non alphanumeric ([a-zA-Z0-9]) characters. 

ieeg_json.SamplingFrequency = 'n/a'; % Sampling frequency (in Hz) of all the iEEG channels 
% in the recording (e.g., 2400). All other channels should have frequency specified as 
% well in the channels.tsv file.

ieeg_json.PowerLineFrequency = 60; % Frequency (in Hz) of the power grid where 
% the iEEG recording was done (i.e. 50 or 60) 

ieeg_json.SoftwareFilters = 'n/a'; % List of temporal software filters applied or ideally 
% key:value pairs of pre-applied filters and their parameter values. (n/a if none).

%% Optional fields:

%ieeg_json.HardwareFilters = 'HighPass:0.16,LowPass:500';% List of hardware (amplifier) filters applied with 
% key:value pairs of filter parameters and their values.
%ieeg_json.HardwareFilters =  sprintf('{"HighPass": {"cutoff":0.16},"LowPass": {"cutoff":500}}');
ieeg_json.HardwareFilters.HighPass.cutoff = 0.16;
ieeg_json.HardwareFilters.LowPass.cutoff = 500;

ieeg_json.Manufacturer = 'Natus'; % Manufacturer of the amplifier system  (e.g. "TDT, blackrock")

ieeg_json.ManufacturersModelName = 'Natus Quantum PK1171'; % Manufacturer's designation of the 
% iEEG amplifier model (e.g. "TDT"). 

ieeg_json.TaskDescription = 'n/a'; % Longer description of the task.

ieeg_json.Instructions = 'n/a'; % Text of the instructions given to participants 
% before the recording. This is especially important in context of resting 
% state and distinguishing between eyes open and eyes closed paradigms. 

% ieeg_json.CogAtlasID = 'n/a'; % URL of the corresponding Cognitive Atlas Task term

% ieeg_json.CogPOID = 'n/a'; %  URL of the corresponding CogPO term

ieeg_json.InstitutionName = 'New York School of Medicine'; %  The name of the institution in charge of 
% the equipment that produced the composite instances.

ieeg_json.InstitutionAddress = '550 1st Avenue, New York, NY 10016, USA'; % The address of the institution 
% in charge of the equipment that produced the composite instances. 

ieeg_json.DeviceSerialNumber = 'n/a'; % The serial number of the equipment that 
% produced the composite instances. A pseudonym can also be used to prevent 
% the equipment from being identifiable, as long as each pseudonym is unique 
% within the dataset.

ieeg_json.ECOGChannelCount = '0'; % Number of iEEG surface channels included in the recording (e.g. 120)

ieeg_json.SEEGChannelCount = '0'; % Number of iEEG depth channels included in the recording (e.g. 8)

ieeg_json.EEGChannelCount = '0'; % Number of EEG channels included in the recording (e.g. 0) 

ieeg_json.EOGChannelCount = '0'; % Number of EOG channels included in the recording (e.g. 1)

ieeg_json.ECGChannelCount = '0'; % Number of ECG channels included in the recording (e.g. 1)

ieeg_json.EMGChannelCount = '0'; % Number of EMG channels included in the recording (e.g. 1)

ieeg_json.MiscChannelCount = '0'; % Number of miscellaneous channels included in 
% the recording (e.g. 1)

ieeg_json.TriggerChannelCount = '0'; % Number of channels for digital (TTL bit level) 
% triggers (e.g. 0) 

ieeg_json.RecordingDuration = '0'; % Length of the recording in seconds (e.g. 3600)

ieeg_json.RecordingType = 'continuous'; % continuous, epoched 

%ieeg_json.EpochLength = 'Inf'; % Duration of individual epochs in seconds (e.g. 1).
% If recording was continuous, set value to Inf.

ieeg_json.SubjectArtefactDescription = 'n/a'; % Freeform description of the observed 
% subject artefact and its possible cause (e.g. door open, nurse walked into room at 2 min, 
% "Vagus Nerve Stimulator", non-removable implant, seizure at 10 min). 
% If this field is left empty, it will be interpreted as absence of artifacts.

ieeg_json.SoftwareVersions = 'Natus Neuroworks'; % Manufacturer's designation of the acquisition software.


%% Specific iEEG fields:


%% Required fields:
ieeg_json.iEEGReference = 'Two contact subdural strips either flipped facing the skull or screwed into the skull'; 
% General description of the reference scheme used and (when applicable) of
% location of the reference electrode in the raw recordings (e.g. "left
% mastoid?, ?bipolar?, ?T01? for electrode with name T01, ?intracranial
% electrode on top of a grid, not included with data?, ?upside down
% electrode?). If different channels have a different reference, this field
% should have a general description and the channel specific reference
% should be defined in the _channels.tsv file.

%% Optional fields:

ieeg_json.ElectrodeManufacturer = 'AdTech'; % Can be used if all electrodes are of the 
% same manufacturer (e.g. AD-TECH, DIXI). If electrodes of different manufacturers 
% are used, please use the corresponding table in the _electrodes.tsv file. 

ieeg_json.ElectrodeManufacturersModelName = 'Standard ad-tech grids (64 8x8 1cm spacing) and strips (1 cm spacing) 4mm contacts, 2.3 mm exposed'; % If different electrode types are used, 
% please use the corresponding table in the _electrodes.tsv file.

ieeg_json.iEEGGround = 'Two contact subdural strips either flipped facing the skull or screwed into the skull'; 
 %Description  of the location of the ground electrode (?placed on right mastoid (M2)?).

ieeg_json.iEEGPlacementScheme = 'n/a'; % General description of the placement 
% of the iEEG electrodes. Left/right/bilateral/depth/surface 
% (e.g. ?left frontal grid and bilateral hippocampal depth? or ?surface strip 
% and STN depth?).

ieeg_json.iEEGElectrodeGroups = 'n/a'; % Field to describe the way electrodes 
% are grouped into strips, grids or depth probes e.g. {'grid1': "10x8 grid 
% on left temporal pole", 'strip2': "1x8 electrode strip on xxx"}.

ieeg_json.ElectricalStimulation = false; % Optional field to specify if electrical stimulation 
% was done during the recording (boolean). Parameters 
% for event-like stimulation should be specified in the _events.tsv file 
% (see example underneath). Continuous parameters that change across ?scans? 
% can be indicated in the the _scans.tsv file.

%ieeg_json.Medication = 'n/a'; %  Optional field to add medication that the patient 
% was on during a recording. 

%% 
json_options.indent = '    '; % this just makes the json file look prettier 
% when opened in a text editor

end