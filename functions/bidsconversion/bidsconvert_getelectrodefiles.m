function [electrode_table, channel_table] = bidsconvert_getelectrodefiles(dataReadDir, hdr, triggerChannel, badChannels, badChannelsDescriptions)

% READ IN relevant data files %%%%%%%%%%%%%%%%%%


%   - electrode coordinate file provided by SoM

% Locate and read the electrode file provided by SOM
% We want the one that matches the T1
ElecFile = dir(fullfile(dataReadDir,'*coor_T1*.txt'));

if length(ElecFile) < 1
    fprintf('[%s] Warning: no coordinate file found!\n', mfilename) 
    elec_labels = hdr.label; chanNames = elec_labels;
    elec_types = hdr.chantype; chanTypes = elec_types;
    chanUnits = hdr.chanunit;
    elec_xyz = nan(length(elec_labels),3);
    elecInx = 1:length(elec_labels);
else
    if length(ElecFile) > 1, fprintf('[%s] Warning: multiple elec files found: using first one', mfilename); end
    fprintf('[%s] Reading electrode information from: %s\n', mfilename, ElecFile(1).name); % if there are multiple coor_T1 files for separate hemispheres, D(1) will always be the full list
    elec_fname = fullfile(ElecFile(1).folder,ElecFile(1).name);
    fid = fopen(elec_fname); E = textscan(fid, '%s%f%f%f%s'); fclose(fid);
    elec_xyz = [E{2} E{3} E{4}]; 
    elec_labels = E{1};
    elec_types = E{5};
    fprintf('[%s] Types of electrodes found: %s\n', mfilename, unique(strcat(elec_types{:})));
    % Replace elec_types with BIDS terminology
    elec_types = replace(elec_types, 'G', 'surface');
    elec_types = replace(elec_types, 'D', 'depth');
    elec_types = replace(elec_types, 'S', 'surface');

    % Match elec names; sort coordinates and electrode types based on matching
    [elecInx, chanNames, chanTypes, chanUnits] = bidsconvert_getchannelspecs(hdr, elec_labels, elec_types);
end

% Generate an  electrodes table with SOM defaults
n = length(elecInx);
[electrode_table] = createBIDS_ieeg_electrodes_tsv_nyuSOM(n);

% Fill with relevant information:
electrode_table.name = elec_labels(elecInx);
electrode_table.x = elec_xyz(elecInx,1);
electrode_table.y = elec_xyz(elecInx,2);
electrode_table.z = elec_xyz(elecInx,3);
electrode_table.type = elec_types(elecInx);

% Generate a channel table with SOM defaults 
[channel_table] = createBIDS_ieeg_channels_tsv_nyuSOM(length(chanNames));
channel_table.name  = chanNames;
channel_table.type  = chanTypes;
if ~isempty(chanUnits), channel_table.units = chanUnits; end
channel_table.sampling_frequency = repmat(hdr.Fs,size(channel_table,1),1);

% Indicate which channel is the trigger channel in channels.tsv
channel_table.type{triggerChannel} = 'trig';

% Indicate channel statuses
%channel_table.status = repmat({'bad'},height(channel_table),1);
%channel_table.status_description = repmat({'bad'},height(channel_table),1);
%channel_table = table2cell(channel_table);
ecogChannels = find(contains(channel_table.type,{'seeg', 'ecog'}));
for ii = 1:length(ecogChannels)
   channel_table.status{ecogChannels(ii)} = 'good';
end
%for ii = 1:length(goodChannels)
%   channel_table.status{goodChannels(ii)} = 'good';
%end
for ii = 1:length(badChannels)
    channel_table.status{badChannels(ii)} = 'bad';
    if exist('badChannelsDescriptions', 'var')
        channel_table.status_description{badChannels(ii)} = badChannelsDescriptions{ii};
    end
end

end