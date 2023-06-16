function [epochs] = ecog_basedetrendEpochs(epochs, t, baselineTime, onsetDelay)
% Applies a baseline detrend to the epoched time courses based on a
% specified time interval within the epoch, and then applies a baseline
% corretion. 
%
% [epochs] = ecog_basedetrendEpochs(epochs, t, baselineTime, onsetDelay)
%
% Input
%   epochs:             3D array containing epoched time series (samples x
%                       epochs x channels)
%   t:                  1D array of length samples indicating time relative 
%                       to stimulus onset (in seconds)
%   baselineTime:       time window over which to compute the baseline
%                       signal (in seconds), e.g [-0.3 0]; 
%                       default: all t smaller than 0
%                       baseline detrend is applied for all t smaller than
%                       upper limit of baselineTime
%   onsetDelay:         1 x epochs array of stimulus onset delay (in seconds)
%
% Output
%   epochs_normalized:  3D array containing normalized epoched time series
%                       (samples x epochs x channels)
% Example
%
% See also ecog_makeEpochs.m, ecog_normalizeEpochs.m


% set default parameter
if ~exist('baselineTime', 'var') || isempty(baselineTime)
    baselineTime = [-inf 0];
elseif numel(baselineTime) < 2
    baselineTime = [-inf baselineTime];
end

% check onset delay
if exist('onsetDelay', 'var') && ~isempty(onsetDelay)
    if numel(onsetDelay)==1
        [epochs] = ecog_basedetrendEpochs(epochs, t-onsetDelay, baselineTime);
    elseif numel(onsetDelay)==size(epochs,2) || numel(onsetDelay)==size(epochs,2:ndims(epochs))
        % reshape in matrix if needed
        epoch_size = size(epochs);
        if numel(onsetDelay)==size(epochs,2:ndims(epochs))
        epochs = reshape(epochs,size(epochs,1),[]);
        else
        epochs = reshape(epochs,size(epochs,1),size(epochs,2),[]);
        end
        % apply baseline detrend for each onset delay
        for idelay = reshape(unique(onsetDelay),1,[])
            groupidx = onsetDelay == idelay;
            [epochs(:,groupidx,:)] = ecog_basedetrendEpochs(epochs(:,groupidx,:), t-idelay, baselineTime);
        end
        % reshape back to origin dims
        epochs = reshape(epochs,epoch_size);        
    else
        error('''onsetDelay'' has invalid size.');
    end
else        % Main function
    % define the reference and baseline window
    ref_range  = t < baselineTime(2);
    base_range = t < baselineTime(2);
    % base_range = t >= baselineTime(1) & t < baselineTime(2);
    base_range(any(isnan(epochs),2:ndims(epochs))) = false;     % omitnan
    
    % reshape in matrix
    epoch_size = size(epochs);
    epochs = reshape(epochs,size(epochs,1),[]);
    
    % detrend
    regref = [[linspace(1,0,sum(ref_range))';zeros(sum(~ref_range),1)], ones(size(t))];  % regressor
    coef = regref(base_range,:) \ epochs(base_range,:);
    epochs = epochs - regref * coef;
    
    % reshape back to origin dims
    epochs = reshape(epochs,epoch_size);
    
    % normalize
    [epochs] = ecog_normalizeEpochs(epochs, t, baselineTime, 'subtractwithintrial');
end

end

