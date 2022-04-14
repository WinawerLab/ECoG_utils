function [electrode_table, channel_table] = bidsconvert_getelectrodefiles(dataReadDir, hdr, triggerChannel, badChannels, badChannelsDescriptions)

% READ IN relevant data files %%%%%%%%%%%%%%%%%%


%   - electrode coordinate file provided by SoM

% Locate and read the electrode file provided by SOM
% We want the one that matches the T1
ElecFile = dir(fullfile(dataReadDir,'*coor_T1*.txt'));

if length(ElecFile) < 1
    warning('No coordinate file found!') 
    elec_labels = hdr.label; chanNames = elec_labels;
    elec_types = hdr.chantype; chanTypes = elec_types;
    
    if length(contains(elec_types, 'unknown')) > 0.9 * length(elec_types)
        % heuristic for chanTypes if no recon available and ft labels nearly all channels as unknowns:
        chanTypes(:) = {'ecog'};
        for ii = 1:length(chanTypes)
            chanNameLengths(ii) = numel(hdr.label{ii});
        end
        ind = contains(hdr.label, 'D'); chanTypes(ind) = {'seeg'};
        ind = contains(hdr.label, {'SG', 'Pleth', 'PR', 'OSAT', 'TRIG'}); chanTypes(ind) = {'other'};
        ind = contains(hdr.label, {'DC'})' & chanNameLengths == 3; chanTypes(ind) = {'other'};
        ind = contains(hdr.label, {'ECG', 'EKG'}); chanTypes(ind) = {'ecg'};
    end

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
    elec_types = replace(elec_types, 'EG', 'surface');
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
electrode_table.hemisphere(electrode_table.x <0) = {'L'};
electrode_table.hemisphere(electrode_table.x >0) = {'R'};

% Generate a channel table with SOM defaults 
[channel_table] = createBIDS_ieeg_channels_tsv_nyuSOM(length(chanNames));
channel_table.name  = chanNames;
channel_table.type  = chanTypes;

% Check the units against the hdr info: use default 'microV' for ieeg
% channels, use the info from the header for other channels; if no
% information in header, put 'unknown'
for ii = 1:length(chanUnits)
    if ~isempty(chanUnits{ii}) && ~strcmp(chanUnits{ii}, 'uV'), channel_table.units{ii} = chanUnits{ii}; end
    if isempty(chanUnits{ii}), channel_table.units{ii} = 'unknown'; end
end
channel_table.sampling_frequency = repmat(hdr.Fs,size(channel_table,1),1);

% Indicate which channel is the trigger channel in channels.tsv
channel_table.type{triggerChannel} = 'trig';

% Indicate which channels below to which group in channels.tsv
INXNames = {'depth', 'grid', 'grid', 'grid', 'grid', 'strip'};
INX = [];
% DEPTH electrodes:
INX{1} = find(contains(lower(channel_table.type), 'seeg'));
% GRID electrodes:
INX{2} = find(contains(lower(channel_table.type), 'ecog') & strncmp('GA', channel_table.name, 2));
INX{3} = find(contains(lower(channel_table.type), 'ecog') & strncmp('GB', channel_table.name, 2));
INX{4} = find(contains(lower(channel_table.type), 'ecog') & strncmp('GC', channel_table.name, 2));  % For future robustness
% if the grid is labeled as 'G', not A or B:
INX{5} = find(contains(lower(channel_table.type), 'ecog') & strncmp('G', channel_table.name,1) & ~strncmp('GA', channel_table.name, 2) & ~strncmp('GB', channel_table.name, 2) & ~strncmp('GC', channel_table.name, 2));
% ALL OTHER SURFACE electrodes
INX{6} = find(contains(lower(channel_table.type), 'ecog') & ~strncmp('G', channel_table.name,1));
% Judge HD grid
girdHDthresh = 64;
INXNames(cellfun(@length,INX)>girdHDthresh & ismember(INXNames,'grid')) = {'HDgrid'};
for ii = 1:length(INX)
    chan_index = INX{ii};
    if ~isempty(chan_index)        
        fprintf('[%s] Assigning electrodes to group: %s...\n',mfilename, INXNames{ii});
        channel_table.group(chan_index,:) = INXNames(ii);
    end
end

% Update electrodes.tsv to have same groups as channel table
elec_inx = ecog_matchChannels(electrode_table.name, channel_table.name);
if length(find(elec_inx>0)) ~= height(electrode_table)
    error('Could not find all electrode names in channel table!')
end
electrode_table.group = channel_table.group(elec_inx);
% For HDgrid electrodes, overwrite the default size
hd_inx = contains(electrode_table.group, 'HDgrid');
electrode_table.size(hd_inx) = 1.0;
electrode_table.manufacturer(hd_inx,:) = {'PMT'};

% Indicate bad channels
for ii = 1:length(badChannels)
    channel_table.status{badChannels(ii)} = 'bad';
    if exist('badChannelsDescriptions', 'var')
        channel_table.status_description{badChannels(ii)} = badChannelsDescriptions{ii};
    end
end

% Also label non-ecog channels as bad
nonEcogChannels = find(~contains(channel_table.type,{'seeg', 'ecog'}));
for ii = 1:length(nonEcogChannels)
   channel_table.status{nonEcogChannels(ii)} = 'bad';
end

% Set channel table type to upper case (bids-requirement)
channel_table.type = upper(channel_table.type);


end