function [epochs] = ecog_normalizeEpochs(epochs, t, baselineTime, baselineType, idx)
% Applies a normalization to the epoched time courses based on a specified
% time interval within the epoch. The normalization can be either a
% conversion to percent signal change (appropriate for broadband) or a
% simple subtraction of the baseline within the epoch itself (appropriate
% for ERPs).
%
% [epochs] = ecog_normalizeEpochs(epochs, t, baselineTime, baselineType, idx)
%
% Input
%   epochs:             3D array containing epoched time series (samples x
%                       epochs x channels)
%   t:                  1D array of length samples indicating time relative 
%                       to stimulus onset (in seconds)
%   baselineTime:       time window over which to compute the baseline
%                       signal (in seconds), e.g [-0.2 0]; 
%                       default: all t smaller than 0
%   baselineType:       string describing the normalization procedure, 
%                       which should be one of the following:
%                           - percentsignalchange (default)
%                           - subtractwithintrial 
%   idx:                vector equal to number of epochs indicating subgroups 
%                       for which to apply normalization separately (e.g.,
%                       individual runs). 
%                       If baselineType issubtractwithintrial, idx is ignored
%                       and nomalization is applied for each epoch.
%                           default: vector of ones (no subgroups).
%
% Output
%   epochs_normalized:  3D array containing normalized epoched time series
%                       (samples x epochs x channels)
% Example
%
% See also ecog_makeEpochs.m

if ~exist('baselineType', 'var') || isempty(baselineType)
    baselineType = 'percentsignalchange';
end

if ~exist('idx', 'var') || isempty(idx)
    idx = ones(size(epochs,2),1); 
end

% define the baseline window
if ~exist('baselineTime', 'var') || isempty(baselineTime)
    base_range = (t < 0);
else
    base_range = (t >= baselineTime(1) & t < baselineTime(2));
end

% determine how many runs we need to normalize for
runs = unique(idx);

% permute dimensions so we can use pointwise division
epochs = permute(epochs, [3 2 1]);

% normalize
for ii = 1:length(runs)
    idx_i = (idx == runs(ii));
    
    switch baselineType
        case 'percentsignalchange'
            %m_base     = squeeze(median(mean(epochs(:,idx_i,base_range),3, 'omitnan'), 2, 'omitnan'));
            m_base     = squeeze(mean(mean(epochs(:,idx_i,base_range),3, 'omitnan'), 2, 'omitnan'));
            epochs(:,idx_i,:) = epochs(:,idx_i,:)./m_base-1;
        case 'subtractwithintrial'
            epochs(:,idx_i,:) = epochs(:,idx_i,:) - mean(epochs(:,idx_i,base_range),3, 'omitnan');
            % debug: compare with Dora's function:
            %epochs(:,idx_i,:) = ecog_baselinesubtract(epochs(:,idx_i,:),base_range);
        otherwise
            error('unknown normalization calculation')
    end
end

% permute back to origin dims
epochs = permute(epochs, [3 2 1]);

end

