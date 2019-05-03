function [dataReadDir, dataWriteDir, stimWriteDir, T1WriteDir] = bidsconvert_getpaths(patientID, RawDataDir, BIDSDataDir, projectName, sub_label, ses_label, ses_labelt1)


if isnumeric(patientID)
    patientID = num2str(patientID);
end

dataReadDir  = fullfile(RawDataDir, patientID);
stimWriteDir = fullfile(BIDSDataDir, projectName, 'stimuli');
dataWriteDir = fullfile(BIDSDataDir, projectName, sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label), 'ieeg');
T1WriteDir   = fullfile(BIDSDataDir, projectName, sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_labelt1), 'anat');

end
