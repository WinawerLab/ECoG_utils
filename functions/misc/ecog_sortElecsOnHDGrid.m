function [chanInx, inxNames] = ecog_sortElecsOnHDGrid(channels,indices)

% Calculates a sorting index so HD electrodes can be plotted according to
% how they are laid out on the 16*8 HD grid used in the BAIR studies
% Index is split in two, for separate plots of bottom and top 

% Electrode positions on the grid
v = 1:16:113;
inx{1} = [v v+1 v+2 v+3 v+4 v+5 v+6 v+7]; % TOP HALF OF HD GRID
inx{2} = inx{1} + 8;                      % BOTTOM HALF OF HD GRID
inxNames = {'top', 'bottom'};

channelNames = channels.name(indices);
chanInx = cell(2,1);

for ii = 1:64
    chanToMatch = sprintf('%s%03d', 'GB', inx{1}(ii)); 
    elecNumber = strcmp(chanToMatch,channelNames);
    if any(elecNumber)
        chanInx{1}(ii) = indices(elecNumber);
    else
        chanInx{1}(ii) = nan; % put empty index for not-connected channel
    end
end

for ii = 65:128
    chanToMatch = sprintf('%s%03d', 'GB', inx{2}(ii-64)); 
    elecNumber = strcmp(chanToMatch,channelNames);
    if any(elecNumber)
        chanInx{2}(ii-64) = indices(elecNumber);
    else
        chanInx{2}(ii-64) = nan; % put empty index for not-connected channel
    end
end

end

% OLD
% gridList = [];
% for ee = 1:128
% 	chanName = sprintf('%s%03d', 'GB', ee);   
%     if strmatch(chanName,gridElecs, 'exact')
%         gridList{ee} = chanName;
%     else
%         gridList{ee} = [chanName '-nodata'];
%     end
% end
% 
% % Electrode positions on the grid
% v = 1:16:113;
% inx{1} = [v v+1 v+2 v+3 v+4 v+5 v+6 v+7]; % TOP HALF OF HD GRID
% inx{2} = inx{1} + 8;                      % BOTTOM HALF OF HD GRID
% inxNames = {'top', 'bottom'};
