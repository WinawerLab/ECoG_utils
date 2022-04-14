function [electrode_table] = createBIDS_ieeg_electrodes_tsv_nyuSOM(n)

% Adapted from BIDS_starter_kit createBIDS_ieeg_electrodes_tsv template
% Intended to be run to get default options for NYU SOM recordings
%
% IGroen, 2018

%% Template Matlab script to create an BIDS compatible electrodes.tsv file
% This example lists all required and optional fields.
% When adding additional metadata please use CamelCase 
%
% DHermes, 2017
% modified RG 201809


%% required columns

name = repmat(' ', n, 1); %{''}; % Name of the electrode
x = zeros(n,1);
y = zeros(n,1);
z = zeros(n,1);
size = ones(n,1)*2.3; % Diameter in mm 
type = repmat(' ', n, 1); % Type of intracranial electrode, one of [?surface?,  ?depth? , ?dbs?] 

%% optional columns

material = repmat({'platinum'}, n, 1); % Recommended. Material of the electrodes

manufacturer = repmat({'AdTech'}, n, 1); % Optional field to specify the electrode manufacturer 
% for each electrode. Can be used if electrodes were manufactured by more than one company.

group = repmat({'n/a'}, n, 1); % Optional field to specify the group that the electrode 
% is a part of. Note that any group specified here should match a group 
% specified in `_channels.tsv` (details of grid size can be added in _ieeg.json 
% the iEEGElectrodeGroups field. )

hemisphere = repmat({'n/a'}, n, 1); % Optional field to specify the hemisphere in which the 
% electrode is placed, one of [L or R] (use capital).

%tissue = repmat('n/a', n, 1); % If a clinician has made an observation about the tissue 
% underneath the electrode  (e.g. epilepsy, tumor, if nothing state n/a)

%% generate table
%electrode_table = table(name,x,y,z,size,type,material, manufacturer);
electrode_table = table(name,x,y,z,size,type,material, manufacturer, group, hemisphere);

end
