function [inx, inxNames] = ecog_sortElecsOnHDGrid(gridElecs)

% Calculates a sorting index so HD electrodes can be plotted according to
% how they are laid out on the 16*8 HD grid used in the BAIR studies
% Index is split in two, for separate plots of bottom and top 

gridList = [];
for ee = 1:128
	chanName = ['GB' num2str(ee)];    
    if strmatch(chanName,gridElecs, 'exact')
        gridList{ee} = chanName;
    else
        gridList{ee} = [chanName '-nodata'];
    end
end

% Plot response per electrode
v = 1:16:113;
inx{1} = [v v+1 v+2 v+3 v+4 v+5 v+6 v+7]; % TOP HALF OF HD GRID
inx{2} = inx{1} + 8;                      % BOTTOM HALF OF HD GRID
inxNames = {'top', 'bottom'};

end