function [channels] = bair_addAnalyzePrfResultsToChannelTable(channels,results, degPerPix)
% Looks up channels by name in analyzePRF results output and adds the
% following parameters to input channels: R2, ecc, ang, rfsize, expt, gain.
%
% [channels] = bair_addAnalyzePrfResultsToChannelTable(channels,results,
% degPerPix)
%
% Input:
%   channels: channels table containing invididual electrode information.
%   results: output file generated using analyzePRF
%   degPerPix: degrees per visual angle (default for BAIR: 16.6/100);
% 
% If 'results' is empty, it will add nans to the channel table.
%
% Note: analyzePRF does not automatically add channel names to the results
% output, so make sure you save those out when fitting the PRFs, as a
% separate field called 'channels' that has a field called 'channels.name'.
%
% IG 2020

if ~exist('degPerPix', 'var') || isempty(degPerPix)
    % The stimulus is 100 pixels (in both height and weight), and this corresponds to
    % 16.6 degrees of visual angle:
    degPerPix= 16.6/100;
end

nChan       = height(channels);
aprf_R2     = nan([nChan  1]); 
aprf_ecc    = nan([nChan  1]); 
aprf_ang    = nan([nChan  1]); 
aprf_rfsize = nan([nChan  1]); 
aprf_expt   = nan([nChan  1]); 
aprf_gain   = nan([nChan  1]); 

if ~isempty(results) && isfield(results, 'channels')
    fprintf('[%s] Adding prf solutions to channel table\n',mfilename);

    for ii = 1:nChan
        % match channels.name to results.channels.name
        chan_idx = find(strcmp(channels.name{ii}, results.channels.name));
        if length(chan_idx) > 1
            warning('[%s] found multiple matches for channel name %s in prf results!\n', mfilename, channels.name{ii});
            continue
        elseif length(chan_idx) < 1
            warning('[%s] did not find any matches for channel name %s in prf results!\n', mfilename, channels.name{ii});
            continue
        end

        % put corresponding values into aprf columns
        aprf_R2(ii)     = results.R2(chan_idx);
        aprf_ecc(ii)    = results.ecc(chan_idx)*degPerPix; 
        aprf_ang(ii)    = results.ang(chan_idx); 
        aprf_rfsize(ii) = results.rfsize(chan_idx)*degPerPix; 
        aprf_expt(ii)   = results.expt(chan_idx); 
        aprf_gain(ii)   = results.gain(chan_idx); 
    end
end

% Add to channel table
channels.aprf_R2     = aprf_R2;
channels.aprf_ecc    = aprf_ecc;
channels.aprf_ang    = aprf_ang;
channels.aprf_rfsize = aprf_rfsize;
channels.aprf_expt   = aprf_expt;
channels.aprf_gain   = aprf_gain;

end