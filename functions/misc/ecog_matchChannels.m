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
           newx = [];
           for kk = 1:length(x)
               if length(chanList{x(kk)}) == length(stringtomatch)
                   newx(kk) = x(kk);
               end            
           end
           if length(newx) == 1
               fprintf('[%s] warning: checked name length and found a single match for %s - continuing\n',mfilename,stringtomatch)
               elInx(ii) = newx;
           else
               fprintf('[%s] warning: checked name length and still found multiple matched for %s, pls check channel names!\n',mfilename,stringtomatch)
               elInx(ii) = x(1);
           end
        else
            elInx(ii) = x(1);
        end     
    end
end

if isempty(elInx)
    fprintf('[%s] could not match this electrode name!\n',mfilename);
    return
end
end
   