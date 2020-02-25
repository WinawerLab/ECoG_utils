function [channels] = bair_addAnalyzePrfResultsToChannelTable(channels,results)


% Check if input channel table matches channel table in aprf results
% add results.ecc, ang, x, y, sd? need to convert from pix (check
% winawerlab sampleData code)
fprintf('[%s] Adding prf solutions to channel table\n',mfilename);

nChan = height(channels);
aprf_R2  = nan([nChan  1]); 
aprf_ecc = nan([nChan  1]); 
aprf_ang = nan([nChan  1]); 

for ii = 1:nChan
    % match channels.name to results.channels.name
    % put corresponding values into aprf columns
    % add to channel table
end

% Add area matches to channel table
channels.aprfR2    = aprf_R2;
channels.aprfeccen = aprf_ecc;
channels.aprfangle = aprf_ang;

end