function [epochs, outlier_idx, max_pows, outlier_thresh] = ecog_selectEpochsStat(epochs, t, stim_on, threshpow, threshdist)
% Epoch selection 
%    
% Remove epochs whose max response power is in <threshdist> significance
% of estimated inverse gaussian distribution in that channel.
% The distribution is estimated from epochs whose max response power is
% <threshpow> times lower than the median response in that channel.

% 20200228 Yuaa: modified from ecog_selectEpochs

%% Set parameters

narginchk(3,5);
if nargin<4,    threshpow  = 20;    end
if nargin<5,    threshdist = 5e-4;  end

%-- Compute time indices
stim_on_idx = t > stim_on(1) & t <= stim_on(2);  

%% 1st step: exclude outliers based on median and variance
%-- Compute max power over time for each epoch
max_pows = shiftdim(nanmax(epochs.^2,[],1),1);

%-- Compute max power over stim_on period across epochs
max_pows_stim_on = shiftdim(nanmax(epochs(stim_on_idx,:,:).^2,[],1),1);

%-- Put all epochs with max higher than outlier_thresh * median to nans
outlier_thresh = threshpow * median(max_pows_stim_on);
outlier_idx = max_pows > outlier_thresh;

% %% 2nd step: estimate chi2 distribution
epochs_conf = epochs; epochs_conf(:,outlier_idx) = nan;
%-- Compute max power over time for each epoch
max_pows_conf = shiftdim(nanmax(epochs_conf.^2,[],1),1);

%-- Fit distribution
%%% fit to InverseGaussian distribution
for el = 1:size(max_pows_conf,2)
    estdist = fitdist(max_pows_conf(:,el),'InverseGaussian');
    outlier_thresh(el) = estdist.icdf(1-threshdist);
end
% %%% fit to empirical distribution
% for el = 1:size(max_pows_conf,2)
%     [estdist,estpoint] = ecdf(max_pows_conf(:,el));
%     outlier_thresh(el) = estpoint(find(estdist>(1-threshdist),1));
% end

outlier_idx = max_pows > outlier_thresh;

%% Output
epochs(:,outlier_idx) = nan;

end
