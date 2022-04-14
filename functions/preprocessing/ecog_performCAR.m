 function [signal_reref, channels_reref, group_indices, group_names] = ecog_performCAR(signal, channels)

% If there is no status column, add one assuming all channels are good.
if ~isfield(summary(channels), 'status') 
    channels.status = repmat({'good'}, [height(channels) 1]);
    fprintf('[%s] No status column found in channels file. Assuming all channels are good.\n', mfilename); 
elseif length(unique(channels.status)) == 1 && contains(unique(channels.status), 'n/a')
    fprintf('[%s] Status column in channels file is n/a. Assuming all channels are good.\n', mfilename); 
    channels.status = repmat({'good'}, [height(channels) 1]);
else
    nBadChannels = length(find(contains(channels.status, 'bad')));
    nGoodChannels = length(find(contains(channels.status, 'good')));
    fprintf('[%s] Found %d good channels and %d bad channels. Including only good channels in CAR.\n', mfilename, nGoodChannels, nBadChannels); 
end

% If there is no group column, use the 'type' column to separate depth and
% surface channels and assume those belong to one group each.
if ~isfield(summary(channels), 'group')
    channels.group = channels.type;
    fprintf('[%s] No group column found in channels file. Using type column to designate groups.\n', mfilename); 
end

% Perform CAR separately for different groups of electrodes as indicated
% by the groups column in the channels tsv.
group_names = unique(channels.group);
group_names = group_names(~contains(lower(group_names), {'n/a', 'unknown', 'other', ...
    'misc', 'ref', 'ecg', 'eog', 'emg', 'trig', 'veog', 'heog', 'audio','pupil', 'eyegaze', 'syslock', 'pd', 'adc', 'dac', 'dbs'})); 
% these are all possible non-ieeg type names allowed by BIDS, in case there
% is no group column and the type column is used to group electrodes.

% Initialize
signal_reref = signal;
channels_reref = channels;
group_indices = cell(size(group_names));

% Loop over electrode groups:
for ii = 1:length(group_names)
    
    fprintf('[%s] Referencing %s electrodes...\n',mfilename, group_names{ii});
    chan_index = find(strcmp(channels.group, group_names{ii}));
    
    if ~isempty(chan_index)
        % Include only good_channels in the common average:
        good_channels = find(contains(channels(chan_index,:).status, 'good'));
        % The mean of the good channels is regressed out of all channels
        % (i.e. including the bad ones)
        temp = ecog_carRegress(signal(chan_index,:), good_channels);    
        signal_reref(chan_index,:) = temp;
        channels_reref.reference(chan_index,:) = {'car'};
    end
    group_indices{ii} = chan_index;
end
end