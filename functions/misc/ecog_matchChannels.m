function elInx = ecog_matchChannels(eltomatch, data)

% elInx = ecog_matchChannels(eltomatch, data)
% 
% matches an input string to the field data.labels (fieldtrip format),
% returns indices into data label


if ~iscell(eltomatch)
    eltomatch = {eltomatch};
end

elInx = zeros(size(eltomatch));

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
        fprintf('[%s] could not locate channel info in data!\n',mfilename);
    end
        
    x = find(strncmp(stringtomatch,chanList,length(stringtomatch)));
    %x = find(strcmp(stringtomatch,chanList));
    if ~isempty(x)
        if length(x) > 1
           fprintf('[%s] warning: found multiple matches for electrode %s!\n',mfilename,stringtomatch);
        end
        elInx(ii) = x(1);
    end
end

if isempty(elInx)
    fprintf('[%s] could not match this electrode name!\n',mfilename);
    return
end
end
   