function [epochs, outlier_idx, max_pows, outlier_thresh] = ecog_selectEpochsStat(epochs, regressed, t, stim_on_time, stim_on_epoch, threshpow, threshdist)
% Epoch selection 
%    
% Remove epochs whose max response power is in <threshdist> significance
% of estimated inverse gaussian distribution in that channel.
% The distribution is estimated from epochs whose max response power is
% <threshpow> times lower than the median response in that channel.

% 20200228 Yuaa: modified from ecog_selectEpochs

%% Set parameters

narginchk(5,7);
if nargin<6,    threshpow  = 20;    end
if nargin<7,    threshdist = 1e-10;  end

%-- Compute time indices
if isempty(stim_on_time)
    stim_on_idx = true(size(t));
else
    stim_on_idx = t > stim_on_time(1) & t <= stim_on_time(2);
end

%% 1st step: exclude outliers based on median of stim_on epochs
%-- Compute max power over time for each epoch
max_pows = shiftdim(nanmax(regressed.^2,[],1),1);

%-- Compute reference max power
max_pows_stim_on    = shiftdim(nanmax(epochs(stim_on_idx,stim_on_epoch,:).^2,[],1),1);

%-- set default for outlier_idx
outlier_idx     = false(size(max_pows_stim_on));
n_outlier       = max(sum(outlier_idx,1));
n_outlier_comp  = -size(outlier_idx,1);

while (n_outlier - n_outlier_comp) > 0
    %-- Compute max power over stim_on period across epochs
    max_pows_conf = max_pows_stim_on;
    max_pows_conf(outlier_idx) = nan;

    %-- Put all epochs with max higher than median
    outlier_thresh     = threshpow * nanmedian(max_pows_conf,1);
    outlier_idx = max_pows_stim_on > outlier_thresh;

    n_outlier_comp  = n_outlier;
    n_outlier       = max(sum(outlier_idx,1));
end

%% 2nd step: exclude outliers based on percentile of stim_on epochs
%-- Compute max power over time for each epoch
max_pows_conf = max_pows_stim_on;       max_pows_conf(outlier_idx) = nan;
%-- Cut upper-tail
outlier_thresh = prctile(max_pows_conf,96,1);
outlier_idx = max_pows_stim_on > outlier_thresh;

%-- Cut lower-tail
outlier_thresh = prctile(max_pows_conf,3,1);
outlier_idx = outlier_idx | (max_pows_stim_on < outlier_thresh);

%% 3rd step: estimate chi2 distribution of stim_on epochs
%-- Compute max power over time for each epoch
max_pows_conf = max_pows_stim_on;       max_pows_conf(outlier_idx) = nan;

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
