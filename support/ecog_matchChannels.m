function elInx = ecog_matchChannels(eltomatch, data)

% matches an input string to the field data.labels (fieldtrip format),
% returns indices

elInx = [];

if ~iscell(eltomatch)
    eltomatch = {eltomatch};
end
for ii = 1:length(eltomatch)
    %stringtomatch = ['EEG ' eltomatch{ii}];
    stringtomatch = eltomatch{ii};
    if isfield(data, 'hdr')
        chanList = data.hdr.label;
    elseif isfield(data, 'label')
        chanList = data.label;
    elseif isfield(data, 'channels')
        chanList = data.channels.name;
    else
        disp('could not locate channel info in data!')
    end
        
    x = find(strncmp(stringtomatch, chanList,length(stringtomatch)));
    if ~isempty(x)
        elInx(ii) = x;
    end
end

if isempty(elInx)
    disp('could not match this electrode name')
    return
end
end
   