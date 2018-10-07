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
    x = find(strncmp(stringtomatch, data.label,length(stringtomatch)));
    if ~isempty(x)
        elInx(ii) = x;
    end
end
if isempty(elInx)
    disp('could not match this electrode name')
    return
end
   