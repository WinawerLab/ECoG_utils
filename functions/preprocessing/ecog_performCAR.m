function [signal_reref, INX, INXNames] = ecog_performCAR(signal, channels)

% If there is no status column, add one assuming all channels are good.
if ~isfield(summary(channels), 'status') 
    channels.status = repmat({'good'}, [height(channels) 1]);
end


% Perform CAR separately for surface and depth electrodes; If there is a HD
% grid or separate clinical grids, perform CAR separately for those

INXNames = {'depth', 'surface grid A', 'surface grid B', 'surface'};
% Define subsets of channels to perform CAR within
INX = [];
% DEPTH electrodes:
INX{1} = find(contains(lower(channels.type), 'seeg'));
% GRID electrodes:
INX{2} = find(contains(lower(channels.type), 'ecog') & strncmp('GA', channels.name, 2));
INX{3} = find(contains(lower(channels.type), 'ecog') & strncmp('GB', channels.name, 2));
% ALL OTHER SURFACE electrodes
INX{4} = find(contains(lower(channels.type), 'ecog') & ~strncmp('GA', channels.name, 2) & ~strncmp('GB', channels.name, 2));
%INX{4} = find(contains(lower(channels.type), 'ecog') & ~contains(channels.name, {'GB', 'GA'}));

signal_reref = signal;
for ii = 1:length(INX)
    chan_index = INX{ii};
    fprintf('[%s] Referencing electrodes of type %s...\n',mfilename, INXNames{ii});

    if ~isempty(chan_index)
        % include only good_channels in the common average:
        good_channels = find(contains(channels(chan_index,:).status, 'good'));
        % the mean of the good channels is regressed out of all channels
        % (i.e. including the bad ones)
        temp = ecog_carRegress(signal(chan_index,:), good_channels);    
        signal_reref(chan_index,:) = temp;
    end
end
end