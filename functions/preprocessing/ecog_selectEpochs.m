function [epochs, outlier_idx, max_epochs] = ecog_selectEpochs(epochs, t, stim_on, thresh)
% Epoch selection 
%    
% Remove epochs whose max response amplitude is thresh x more
% or less than the median response in that channel

% Compute time indices
stim_on_idx = t > stim_on(1) & t <= stim_on(2);  

% Compute max over time for each epoch
max_epochs = squeeze(max(epochs,[],1));

% Compute max over stim_on period across epochs
max_epochs_stim_on = squeeze(max(epochs(stim_on_idx,:,:),[],1));

% Put all epochs with max higher than outlier_thresh * median to nans
outlier_thresh = thresh * median(max_epochs_stim_on);
outlier_idx = max_epochs > outlier_thresh;             
epochs(:,outlier_idx) = nan;

end
