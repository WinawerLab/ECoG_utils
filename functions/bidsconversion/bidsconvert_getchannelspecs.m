function [elecInx, chanNames, chanTypes, chanUnits] = bidsconvert_getChannelSpecs(hdr, elecList, elecTypes)

% This function compares channel info from the data (obtained from the file
% header hdr) with electrode info (elecList, obtained from electrode
% coordinate file provided by SOM) to generate the following BIDS
% compatible outputs:
%
% - elecInx  : for each matching channel in the datafile the corresponding
%              index into the electrode file to be used for electrodes.tsv 
% - chanNames: list of shortened names to be used for channels.tsv 
% - chanTypes: list of types to be used for channels.tsv
% - chanUnits: reads the chanunits from the header for channels.tsv

if nargin < 3 || isempty(elecTypes)
    elecTypes = cell(size(elecList));
end

elecInx = nan(length(hdr.label),1); % this will be an index into elecList
chanNames = hdr.label;
chanTypes = hdr.chantype;

if ~isequal(unique(hdr.chanunit), {'unknown'})
    chanUnits = hdr.chanunit;
else
    chanUnits = [];
end

% Match channel names from data header with names from electrode file

if max(contains(hdr.label, 'REF')) == 1
    % this is an older SOM dataset in which channels have not been renamed
    for ii = 1:length(hdr.label)
        C = strsplit(hdr.label{ii}, {' ', '_', '-'});
        % get rid of -REF field
        C = C(1:end-1);
        switch C{1}
            case 'EEG'
                chanName = strcat(C{2:end});
                %inx = find(contains(elecList, chanName));
                inx = find(strncmp(elecList, chanName, length(chanName)));
                if length(inx) == 1
                    elecInx(ii) = inx;
                elseif length(inx) < 1
                    fprintf('[%s] Warning: could not find coordinates for channel %s with electrode information! \n', mfilename, chanName);
                    chanNames{ii} = chanName;
                    chanTypes{ii} = 'unknown';
                elseif length(inx) > 1
                    error('[%s] Multiple coordinates found for channel %s! \n', mfilename, chanName);
                end
                chanNames{ii} = chanName;
            case 'ECG'   
                chanName = strcat(C{2:end});
                chanNames{ii} = chanName;
                chanTypes{ii} = 'ecg';
            otherwise
                chanNames{ii} = strcat(C{:});
                chanTypes{ii} = 'other'; 
        end
    end
else
	% this is a newer SOM dataset in which channels been renamed to no
	% longer have 'REF' or dahes, but the numbering style is different in
	% the elecfile and the hdr so we still need to rename. In this case
	% we'll rename the elecList (remove the 0s):
    elecList_renamed = elecList;
    for ii = 1:length(elecList)
        zeroIdx = strfind(elecList{ii}, '0'); % find the location of 0
        lastIdx = length(elecList{ii}); % check if location of 0 is not at the end
        if ~isempty(zeroIdx) & zeroIdx ~= lastIdx & isempty(str2num(elecList{ii}(zeroIdx-1)))
            keepIdx = setdiff(1:lastIdx, zeroIdx);
            elecList_renamed{ii} = elecList{ii}(keepIdx);
        end
    end
    % Now, loop across all the channels to get their elecInx
    for ii = 1:length(hdr.label)
        chanName = hdr.label{ii};
        if contains(chanName, {'EKG', 'ECG'})
            % this is an ECG channel
            chanNames{ii} = chanName;
            chanTypes{ii} = 'ecg';
        elseif contains(chanName, {'DC', 'TRIG', 'Pleth', 'PR', 'OSAT'})
            % this is a DC channel
            chanNames{ii} = chanName;
            chanTypes{ii} = 'other';
        else
            inx = find(strcmp(chanName,elecList_renamed));
            if length(inx) == 1
                elecInx(ii) = inx;
                chanNames{ii} = elecList{inx};
            elseif length(inx) < 1
                fprintf('[%s] Warning: could not find coordinates for channel %s! \n', mfilename, chanName);
                chanNames{ii} = chanName;
                chanTypes{ii} = 'other';
            elseif length(inx) > 1
                error('[%s] Multiple coordinates found for channel %s! \n', mfilename, chanName);
            end
        end
    end
end

for ii = 1:length(elecInx)
    if ~isnan(elecInx(ii))
        switch elecTypes{elecInx(ii)}
            case 'depth'
                chanTypes{ii} = 'seeg';
            case 'surface'
                chanTypes{ii} = 'ecog';
            otherwise
                chanTypes{ii} = 'other';
        end
    end
end

elecInx = elecInx(~isnan(elecInx));
fprintf('[%s] Found coordinates for %d out of %d electrodes and %d channels \n', mfilename, length(elecInx), length(elecList), length(hdr.label));



end