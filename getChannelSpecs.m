function [elecInx, chanNames, chanTypes, chanUnits] = getChannelSpecs(hdr, elecList)

% This function compares channel info from the data (chanList, obtained
% from the file header) with electrode info (elecList, obtained from
% electrode coordinate file provided by SOM) to generate the following BIDS
% compatible outputs:
% - elecInx: for each matching channel in the datafile the corresponding
% index into the electrode file (to be used for electrodes.tsv)
% - chanNames: list of shortened names to be used for channels.tsv 
% - chanTypes: list of types to be used for channels.tsv

elecInx = nan(length(hdr.label),1); % this will be an index into elecList
chanNames = hdr.label;
chanTypes = hdr.chantype;
chanUnits = hdr.chanunit;

% Match channel names from data header with names from electrode file
for ii = 1:length(hdr.label)
    C = strsplit(hdr.label{ii}, {' ', '_', '-'});
    % get rid of -REF field
    C = C(1:end-1);
    switch C{1}
        case 'EEG'
            chanName = strcat(C{2:end});
            inx = find(contains(elecList, chanName));
            if length(inx) == 1
                elecInx(ii) = inx;
            elseif length(inx) < 1
                fprintf('Warning: could not match SEEG channel %s! \n', chanName);
            elseif length(inx) > 1
                fprintf('Warning: multiple matches with SEEG channel %s! \n', chanName);
            end
            chanNames{ii} = chanName;
            chanTypes{ii} = 'seeg';
        case 'ECG'   
            chanName = strcat(C{2:end});
            chanNames{ii} = chanName;
            chanTypes{ii} = 'ecg';
        otherwise
            chanNames{ii} = strcat(C{:});
            chanTypes{ii} = 'other'; 
    end
end

elecInx = elecInx(~isnan(elecInx));
fprintf('Matched %d out of %d channels and %d electrodes \n', length(elecInx), length(hdr.label), length(elecList));



end