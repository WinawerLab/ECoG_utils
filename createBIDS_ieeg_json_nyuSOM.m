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

ieeg_json.TaskName = ''; % Name of the task (for resting state use the rest
% prefix). No two tasks should have the same name. Task label is derived 
% from this field by removing all non alphanumeric ([a-zA-Z0-9]) characters. 

ieeg_json.SamplingFrequency = ''; % Sampling frequency (in Hz) of all the iEEG channels 
% in the recording (e.g., 2400). All other channels should have frequency specified as 
% well in the channels.tsv file.

ieeg_json.PowerLineFrequency = '60'; % Frequency (in Hz) of the power grid where 
% the iEEG recording was done (i.e. 50 or 60) 

%% Optional fields:

ieeg_json.Manufacturer = ''; % Manufacturer of the amplifier system  (e.g. "TDT, blackrock")

ieeg_json.ManufacturersModelName = ''; % Manufacturer's designation of the 
% iEEG amplifier model (e.g. "TDT"). 

ieeg_json.TaskDescription = ''; % Longer description of the task.

ieeg_json.Instructions = ''; % Text of the instructions given to participants 
% before the recording. This is especially important in context of resting 
% state and distinguishing between eyes open and eyes closed paradigms. 

ieeg_json.CogAtlasID = ''; % URL of the corresponding Cognitive Atlas Task term

ieeg_json.CogPOID = ''; %  URL of the corresponding CogPO term

ieeg_json.InstitutionName = 'New York School of Medicine'; %  The name of the institution in charge of 
% the equipment that produced the composite instances.

ieeg_json.InstitutionAddress = '550 1st Avenue, New York, NY 10016, USA'; % The address of the institution 
% in charge of the equipment that produced the composite instances. 

ieeg_json.DeviceSerialNumber = 'unknown'; % The serial number of the equipment that 
% produced the composite instances. A pseudonym can also be used to prevent 
% the equipment from being identifiable, as long as each pseudonym is unique 
% within the dataset.

ieeg_json.ECOGChannelCount = ''; % Number of iEEG surface channels included in the recording (e.g. 120)

ieeg_json.SEEGChannelCount = ''; % Number of iEEG depth channels included in the recording (e.g. 8)

ieeg_json.EEGChannelCount = ''; % Number of EEG channels included in the recording (e.g. 0) 

ieeg_json.EOGChannelCount = ''; % Number of EOG channels included in the recording (e.g. 1)

ieeg_json.ECGChannelCount = ''; % Number of ECG channels included in the recording (e.g. 1)

ieeg_json.EMGChannelCount = ''; % Number of EMG channels included in the recording (e.g. 1)

ieeg_json.MiscChannelCount = ''; % Number of miscellaneous channels included in 
% the recording (e.g. 1)

ieeg_json.TriggerChannelCount = ''; % Number of channels for digital (TTL bit level) 
% triggers (e.g. 0) 

ieeg_json.RecordingDuration = ''; % Length of the recording in seconds (e.g. 3600)

ieeg_json.RecordingType = 'continuous'; % continuous, epoched 

ieeg_json.EpochLength = 'Inf'; % Duration of individual epochs in seconds (e.g. 1).
% If recording was continuous, set value to Inf.

ieeg_json.SubjectArtefactDescription = ''; % Freeform description of the observed 
% subject artefact and its possible cause (e.g. door open, nurse walked into room at 2 min, 
% "Vagus Nerve Stimulator", non-removable implant, seizure at 10 min). 
% If this field is left empty, it will be interpreted as absence of artifacts.

ieeg_json.SoftwareVersions = ''; % Manufacturer's designation of the acquisition software.


%% Specific iEEG fields:


%% Required fields:
ieeg_json.iEEGReference = ''; % General description of the reference scheme used and (when applicable) of
% location of the reference electrode in the raw recordings (e.g. "left
% mastoid?, ?bipolar?, ?T01? for electrode with name T01, ?intracranial
% electrode on top of a grid, not included with data?, ?upside down
% electrode?). If different channels have a different reference, this field
% should have a general description and the channel specific reference
% should be defined in the _channels.tsv file.


%% Optional fields:

ieeg_json.HardWareFilters = ''; % List of hardware (amplifier) filters applied with 
% key:value pairs of filter parameters and their values.

ieeg_json.SoftWareFilters = ''; % List of temporal software filters applied or ideally 
% key:value pairs of pre-applied filters and their parameter values. (n/a if none).

ieeg_json.ElectrodeManufacturer = ''; % Can be used if all electrodes are of the 
% same manufacturer (e.g. AD-TECH, DIXI). If electrodes of different manufacturers 
% are used, please use the corresponding table in the _electrodes.tsv file. 

ieeg_json.ElectrodeManufacturersModelName = ''; % If different electrode types are used, 
% please use the corresponding table in the _electrodes.tsv file.

ieeg_json.iEEGGround = ''; %Description  of the location of the ground electrode 
% (?placed on right mastoid (M2)?).

ieeg_json.iEEGPlacementScheme = ''; % General description of the placement 
% of the iEEG electrodes. Left/right/bilateral/depth/surface 
% (e.g. ?left frontal grid and bilateral hippocampal depth? or ?surface strip 
% and STN depth?).

ieeg_json.iEEGElectrodeGroups = ''; % Field to describe the way electrodes 
% are grouped into strips, grids or depth probes e.g. {'grid1': "10x8 grid 
% on left temporal pole", 'strip2': "1x8 electrode strip on xxx"}.

ieeg_json.Stimulation = ''; % Optional field to specify if electrical stimulation 
% was done during the recording (options are 1 for yes, 0 for no). Parameters 
% for event-like stimulation should be specified in the _events.tsv file 
% (see example underneath). Continuous parameters that change across ?scans? 
% can be indicated in the the _scans.tsv file.

ieeg_json.Medication = ''; %  Optional field to add medication that the patient 
% was on during a recording. 

%% write
% jsonSaveDir = fileparts(ieeg_json_name);
% if ~isdir(jsonSaveDir)
%     fprintf('Warning: directory to save json file does not exist, create: %s \n',jsonSaveDir)
% end

json_options.indent = '    '; % this just makes the json file look prettier 
% when opened in a text editor

%jsonwrite(ieeg_json_name,ieeg_json,json_options)

end