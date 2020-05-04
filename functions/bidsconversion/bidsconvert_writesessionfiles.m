function bidsconvert_writesessionfiles(dataReadDir, dataWriteDir, T1WriteDir, ...
    sub_label, ses_label, acq_label, ses_labelt1, ...
    electrode_table, dataFileNames, runTimes)


%% Create SESSION-SPECIFIC files %%%%%%%%%%%%%%%%%%
%
%   - T1w.nii.gz
%   - coordsystem.json
%   - electrodes.tsv
%   - scans.tsv
%
% Should be run AFTER generating the run specific files (so we can generate
% the scans.tsv files based on the existing run data files)

%% Create T1w file

% Locate T1 provided by SoM
T1File = dir(fullfile(dataReadDir, 'T1.nii.gz'));
if length(T1File) > 1, fprintf('[%s] Warning: multiple T1s found: using first one\n', mfilename); end
writeT1json = 1;

% Copy file: rename and save in appropriate folder
if length(T1File) < 1
    disp('Warning: no T1 found!') 
    writeT1json = 0;
else
    T1_name = sprintf('sub-%s_ses-%s_acq-%s_run-01_T1w.nii.gz', sub_label, ses_labelt1, acq_label);
    T1_fname = fullfile(T1WriteDir, T1_name);
    if exist(T1_fname, 'file')
        fprintf('[%s] T1 already exists in %s, not overwriting.\n', mfilename, T1_fname);
    else
        fprintf('[%s] Writing %s \n', mfilename, T1_fname);
        copyfile(fullfile(dataReadDir, 'T1.nii.gz'), T1_fname);
        
        % Generate scans table
        filename = string(fullfile('anat', T1_name));
        acq_time = "n/a"; % We don't know when the T1 was done
        T1_scans_table = table(filename, acq_time);
        % Generate output name (hacky but we need to go one level up from anat)
        tmp = strsplit(fileparts(T1WriteDir), '/');
        T1sessionDir = strjoin(tmp,'/');
        T1_scans_name = fullfile(T1sessionDir, sprintf('sub-%s_ses-%s_scans.tsv', sub_label, ses_labelt1));
        writetable(T1_scans_table,T1_scans_name,'FileType','text','Delimiter','\t');

    end
end

%% Create coordsystem.json file
if writeT1json
    % Generate output name
    coord_json_name = fullfile(dataWriteDir, sprintf('sub-%s_ses-%s_acq-%s_coordsystem.json', sub_label, ses_label, acq_label));

    % Get default values
    [coordsystem_json, json_options] = createBIDS_ieeg_coordsystem_json_nyuSOM();

    % Update values with specs for this subject
    coordsystem_json.IntendedFor = fullfile(sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_labelt1), 'anat', ...
        T1_name); % this path must be specified relative to the project folder

    % Write coordsystem.json file
    fprintf('[%s] Writing %s\n', mfilename, coord_json_name);
    jsonwrite(coord_json_name,coordsystem_json,json_options);
end

% Note on coordystem json file: SOM also provides electrode locations in
% MNI space, as well as a list of matched anatomical regions. Do we want to
% include this information in our BIDS dataset as well? The MNI coordinates
% could be included as another coordsystem file (and associated
% electrodes.tsv file) which can be distinguished from the T1 coordinates
% by adding a [_space-<label>] to the filenames (see BIDS iEEG section
% 3.4.1). SOM does not provide an MNI aligned T1 to go along with this file

%% Create electrodes.tsv file

% Generate output name
electrodes_tsv_name = fullfile(dataWriteDir, sprintf('sub-%s_ses-%s_acq-%s_electrodes.tsv', sub_label, ses_label, acq_label));

% Write tsv file
fprintf('[%s] Writing %s\n', mfilename, electrodes_tsv_name);
writetable(electrode_table,electrodes_tsv_name,'FileType','text','Delimiter','\t');

%% Create scans.tsv files

% ieeg

% Generate scans table
acq_time = cell(length(runTimes),1);
for ii = 1:length(runTimes)
    t = runTimes{ii};
    t = datestr(datenum(t, 'yyyymmddTHHMMSS'), 'yyyy-mm-ddTHH:MM:SS');
    % shift date
    t(1:4) = '1900';
    acq_time{ii} = t;
end
filename = dataFileNames;
ieeg_scans_table = table(filename, acq_time);
% Generate output name (hacky but we need to go one level up from ieeg)
tmp = strsplit(fileparts(dataWriteDir), '/');
ieegsessionDir = strjoin(tmp,'/');
ieeg_scans_name = fullfile(ieegsessionDir, sprintf('sub-%s_ses-%s_scans.tsv', sub_label, ses_label));
writetable(ieeg_scans_table,ieeg_scans_name,'FileType','text','Delimiter','\t');

end