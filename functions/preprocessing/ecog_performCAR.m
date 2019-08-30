function [signal_reref, channels_reref, group_indices, group_names] = ecog_performCAR(signal, channels)

% If there is no status column, add one assuming all channels are good.
if ~isfield(summary(channels), 'status') 
    channels.status = repmat({'good'}, [height(channels) 1]);
    fprintf('[%s] No status column found in channels file. Using all channels for CAR.\n', mfilename); 
else
    nBadChannels = length(find(contains(channels.status, 'bad')));
    nGoodChannels = length(find(contains(channels.status, 'good')));
    fprintf('[%s] Found %d good channels and %d bad channels. Including only good channels in CAR.\n', mfilename, nGoodChannels, nBadChannels); 
end


% Perform CAR separately for surface and depth electrodes; If there is a HD
% grid or separate clinical grids, perform CAR separately for those

% INXNames = {'depth', 'surfacegridA', 'surfacegridB', 'surfacestrips'};
% % Define subsets of channels to perform CAR within
% INX = [];
% % DEPTH electrodes:
% INX{1} = find(contains(lower(channels.type), 'seeg'));
% % GRID electrodes:
% INX{2} = find(contains(lower(channels.type), 'ecog') & strncmp('GA', channels.name, 2));
% INX{3} = find(contains(lower(channels.type), 'ecog') & strncmp('GB', channels.name, 2));
% % ALL OTHER SURFACE electrodes
% INX{4} = find(contains(lower(channels.type), 'ecog') & ~strncmp('GA', channels.name, 2) & ~strncmp('GB', channels.name, 2));
% %INX{4} = find(contains(lower(channels.type), 'ecog') & ~contains(channels.name, {'GB', 'GA'}));

group_names = unique(channels.group);
group_names = group_names(~contains(group_names, {'n/a', 'unknown'}));

signal_reref = signal;
channels_reref = channels;
group_indices = cell(size(group_names));

for ii = 1:length(group_names)
    
    fprintf('[%s] Referencing electrodes of group %s...\n',mfilename, group_names{ii});
    chan_index = find(contains(channels.group, group_names{ii}));
    
    if ~isempty(chan_index)
        % include only good_channels in the common average:
        good_channels = find(contains(channels(chan_index,:).status, 'good'));
        % the mean of the good channels is regressed out of all channels
        % (i.e. including the bad ones)
        temp = ecog_carRegress(signal(chan_index,:), good_channels);    
        signal_reref(chan_index,:) = temp;
        channels_reref.reference(chan_index,:) = {'car'};
    end
    group_indices{ii} = chan_index;
end
end