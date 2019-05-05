function [dataReadDir, dataWriteDir, stimWriteDir, T1WriteDir, preprocDir] = ...
    bidsconvert_getpaths(patientID, RawDataDir, BIDSDataDir, projectName, sub_label, ses_label, ses_labelt1)


if isnumeric(patientID)
    patientID = num2str(patientID);
end

dataReadDir    = fullfile(RawDataDir, patientID);
stimWriteDir   = fullfile(BIDSDataDir, projectName, 'stimuli');
dataWriteDir   = fullfile(BIDSDataDir, projectName, sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label), 'ieeg');
T1WriteDir     = fullfile(BIDSDataDir, projectName, sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_labelt1), 'anat');
preprocDir     = fullfile(BIDSDataDir, projectName, 'derivatives', 'preprocessed', sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label));

% Check whether we have write directories
if ~exist(dataWriteDir, 'dir'); mkdir(dataWriteDir);end
if ~exist(stimWriteDir, 'dir'); mkdir(stimWriteDir);end
if ~exist(T1WriteDir, 'dir'); mkdir(T1WriteDir);end
if ~exist(preprocDir, 'dir'); mkdir(preprocDir);end
if ~exist(fullfile(preprocDir, 'figures', 'bidsconversion'), 'dir'); mkdir(fullfile(preprocDir, 'figures', 'bidsconversion')); end

end
