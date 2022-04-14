function [elecInx, chanNames, chanTypes, chanUnits] = bidsconvert_getchannelspecs(hdr, elecList, elecTypes)

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

if ~isequal(unique(hdr.chanunit), {'other'})
    chanUnits = hdr.chanunit;
else
    chanUnits = [];
end

% Match channel names from data header with names from electrode file

if max(contains(hdr.label, 'REF')) == 1
    % this is an older SOM dataset in which channels have not been renamed
    for ii = 1:length(hdr.label)
        C = strsplit(hdr.label{ii}, {' ', '_', '-'});
        if size(C{end},2) == 3
            % get rid of -REF field
            C = C(1:end-1);
        else
            C{end} = C{end}(1:end-3);
            C = {[C{:}]};
        end
        if ~isempty(C)
            if contains(C{1}, 'EEG')
                if max(size(C)) == 1
                    chanName = C{1}(4:end);
                else
                    chanName = strcat(C{2:end});
                end
                %inx = find(contains(elecList, chanName));
                inx = find(strncmp(elecList, chanName, length(chanName)));
                if length(inx) == 1
                    elecInx(ii) = inx;
                elseif length(inx) < 1                   
                    warning('[%s] Could not find coordinates for channel %s in electrode file!', mfilename, chanName);
                    chanNames{ii} = chanName;
                    chanTypes{ii} = 'other';
                elseif length(inx) > 1
                    error('[%s] Multiple coordinates found for channel %s!', mfilename, chanName);
                end
                chanNames{ii} = chanName;
            elseif contains (C{1}, 'ECG')
                if max(size(C)) == 1
                    chanName = C{1}(4:end);
                else
                    chanName = strcat(C{2:end});
                end
                chanNames{ii} = chanName;
                chanTypes{ii} = 'ecg';
            else
                chanNames{ii} = strcat(C{:});
                chanTypes{ii} = 'other'; 
            end
        else
            chanName{ii} = hdr.label{ii};
            chanTypes{ii} = 'other'; 
        end        
    end
    if size(chanTypes,2) > size(chanTypes,1)
        chanTypes = chanTypes';
    end
else
	% this is a newer SOM dataset in which channels been renamed to no
	% longer have 'REF' or dashes, but the numbering style is different in
	% the elecfile and the hdr so we still need to rename. In this case
	% we'll rename the elecList (remove the 0s):
    elecList_renamed = elecList;
    for ii = 1:length(elecList)
        zeroIdx = strfind(elecList{ii}, '0'); % find the location of 0       
        if ~isempty(zeroIdx) 
            keepIdx = setdiff(1:length(elecList{ii}), zeroIdx);
            % double check if we indeed want to remove these zeros or not
            for jj = 1:length(zeroIdx)
                % if location of 0 is at the end, do not remove
                if zeroIdx(jj) == length(elecList{ii}), keepIdx = [keepIdx zeroIdx(jj)]; end
                % if position before 0 is numeric and not zero, do not remove
                tmp = str2num(elecList{ii}(zeroIdx(jj)-1));
                if isnumeric(tmp) & tmp ~= 0, keepIdx = [keepIdx zeroIdx(jj)]; end
            end
            elecList_renamed{ii} = elecList{ii}(unique(keepIdx));
        end
    end
    % Now, loop across all the channels to get their elecInx
    for ii = 1:length(hdr.label)
        chanName = hdr.label{ii};
        if contains(chanName, {'EKG', 'ECG'})
            % this is an ECG channel
            chanNames{ii} = chanName;
            chanTypes{ii} = 'ecg';
        elseif strcmp(chanName, 'DC') | strcmp(chanName, 'TRIG') | strcmp(chanName, 'Pleth') | strcmp(chanName, 'PR') | strcmp(chanName, 'OSAT')
            % this is a DC or other status channel
            chanNames{ii} = chanName;
            chanTypes{ii} = 'other';
        else
            inx = find(strcmp(chanName,elecList_renamed));
            if length(inx) == 1
                elecInx(ii) = inx;
                chanNames{ii} = elecList{inx};
            elseif length(inx) < 1
                warning('[%s]  could not find coordinates for channel %s!', mfilename, chanName);
                chanNames{ii} = chanName;
                chanTypes{ii} = 'other';
            elseif length(inx) > 1
                error('[%s] Multiple coordinates found for channel %s!', mfilename, chanName);
            end
        end
    end
end

for ii = 1:length(elecInx)
    if ~isnan(elecInx(ii))
        if contains(elecTypes{elecInx(ii)}, 'depth')
            chanTypes{ii} = 'seeg';
        elseif contains(elecTypes{elecInx(ii)}, 'surface')
            chanTypes{ii} = 'ecog';
        else
            chanTypes{ii} = 'other';
        end
    end
end

elecInx = elecInx(~isnan(elecInx));
fprintf('[%s] Found coordinates for %d out of %d electrodes and %d channels \n', mfilename, length(elecInx), length(elecList), length(hdr.label));



end