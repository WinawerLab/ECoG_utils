function bidsconvert_writesessionfiles(dataReadDir, dataWriteDir, T1WriteDir, sub_label, ses_label, ses_labelt1, electrode_table)


%% Create SESSION-SPECIFIC files %%%%%%%%%%%%%%%%%%

%   - T1w.nii.gz
%   - coordsystem.json
%   - electrodes.tsv

%% Create T1w file

% Locate T1 provided by SoM
T1File = dir(fullfile(dataReadDir, 'T1.nii.gz'));
if length(T1File) > 1, fprintf('[%s] Warning: multiple T1s found: using first one\n', mfilename); end

% Copy file: rename and save in appropriate folder
if length(T1File) < 1
    disp('Warning: no T1 found!') 
else
    T1_name = fullfile(T1WriteDir, sprintf('sub-%s_ses-%s_T1w.nii.gz', sub_label, ses_labelt1));
    if exist(T1_name, 'file')
        fprintf('[%s] T1 already exists in %s, not overwriting.\n', mfilename, T1_name);
    else
        fprintf('[%s] Writing %s \n', mfilename, T1_name);
        copyfile(fullfile(dataReadDir, 'T1.nii.gz'), T1_name);
    end
end

% Note on pial file: surface recons should go in derivatives folder (see
% BIDS spec doc). If we want to include this in our data do we use our own
% freesurfer recon (which we use for electrode plotting) or the one
% provided by SOM (but why is that a .mat file)? If so, we should probably
% also provide an electrode locations file based on matched nodes
% (see run_E2N and electrode_to_nearest_node scripts in ECoG_utils)

%% Create coordsystem.json file

% Generate output name
coord_json_name = fullfile(dataWriteDir, sprintf('sub-%s_ses-%s_coordsystem.json', sub_label, ses_label));

% Get default values
[coordsystem_json, json_options] = createBIDS_ieeg_coordsystem_json_nyuSOM();

% Update values with specs for this subject
coordsystem_json.IntendedFor = fullfile(sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_labelt1), 'anat', ...
    sprintf('sub-%s_ses-%s_T1w.nii.gz', sub_label, ses_labelt1)); % this path must be specified relative to the project folder

% Write coordsystem.json file
fprintf('[%s] Writing %s\n', mfilename, coord_json_name);
jsonwrite(coord_json_name,coordsystem_json,json_options);

% Note on coordystem json file: SOM also provides electrode locations in
% MNI space, as well as a list of matched anatomical regions. Do we want to
% include this information in our BIDS dataset as well? The MNI coordinates
% could be included as another coordsystem file (and associated
% electrodes.tsv file) which can be distinguished from the T1 coordinates
% by adding a [_space-<label>] to the filenames (see BIDS iEEG section
% 3.4.1). SOM does not provide an MNI aligned T1 to go along with this file

%% Create electrodes.tsv file

% Generate output name
electrodes_tsv_name = fullfile(dataWriteDir, sprintf('sub-%s_ses-%s_electrodes.tsv', sub_label, ses_label));

% Write tsv file
fprintf('[%s] Writing %s\n', mfilename, electrodes_tsv_name);
writetable(electrode_table,electrodes_tsv_name,'FileType','text','Delimiter','\t');

end