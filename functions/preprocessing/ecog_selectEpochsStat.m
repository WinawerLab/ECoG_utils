function [epochs, outlier_idx, max_pows, outlier_thresh] = ecog_selectEpochsStat(epochs, regressed, t, opts)
% Epoch selection 
%    
% Remove epochs whose max response power is in <opts.epoch_outlier_dist>
% significance of estimated inverse gaussian distribution in that channel.
% The distribution is estimated from epochs whose max response power is
% <opts.epoch_outlier_pow> times lower than the median response in that channel.

% 20200228 Yuasa: modified from ecog_selectEpochs
% 20220215 Yuasa: enable to set time range to check epochs

%% Set parameters

narginchk(3,4);
SetDefault('opts.stim_on_time',[]);
SetDefault('opts.stim_on_epoch',true(size(epochs,2)));
SetDefault('opts.epoch_time',[]);
SetDefault('opts.epoch_jump_thresh',250);
SetDefault('opts.epoch_outlier_pow',20);
SetDefault('opts.epoch_outlier_dist',1e-5);

%-- Set parameters
stim_on_time    = opts.stim_on_time;
stim_on_epoch   = opts.stim_on_epoch;
epoch_time      = opts.epoch_time;
jumpthresh      = opts.epoch_jump_thresh;
threshdist      = opts.epoch_dist_thresh;
threshout       = opts.epoch_outlier_thresh;

%-- Compute time indices
if isempty(stim_on_time)
    stim_on_idx = true(size(t));
else
    stim_on_idx = t > stim_on_time(1) & t <= stim_on_time(2);
end
if isempty(epoch_time)
    epoch_idx = true(size(t));
else
    epoch_idx = t > epoch_time(1) & t <= epoch_time(2);
end

%% 1st step: Remove epochs in which there is a sudden jump in voltage from one timepoint to the next.

diff_epochs = abs(diff(epochs(epoch_idx,:,:), 1, 1));
outlier_idx1 = shiftdim(any(diff_epochs > jumpthresh, 1),1);

%-- Set nan for outlier
epochs(:,outlier_idx1)    = nan;
regressed(:,outlier_idx1) = nan;

%% 2nd step: exclude outliers based on distribution of maximum power of regressed signal
%% 2-1 step: temporaly exclude outliers from stim_on distribution based on median of stim_on epochs
%-- Compute max power over time for each epoch
max_pows = shiftdim(nanmax(regressed(epoch_idx,:,:).^2,[],1),1);

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
    outlier_thresh     = threshdist * nanmedian(max_pows_conf,1);
    outlier_idx = max_pows_stim_on > outlier_thresh;

    n_outlier_comp  = n_outlier;
    n_outlier       = max(sum(outlier_idx,1));
end

%% 2-2 step: exclude outliers from stim_on distribution based on percentile of stim_on epochs
%-- Compute max power over time for each epoch
max_pows_conf = max_pows_stim_on;       max_pows_conf(outlier_idx) = nan;
%-- Cut upper-tail
outlier_thresh = prctile(max_pows_conf,96,1);
outlier_idx = max_pows_stim_on > outlier_thresh;

%-- Cut lower-tail
outlier_thresh = prctile(max_pows_conf,3,1);
outlier_idx = outlier_idx | (max_pows_stim_on < outlier_thresh);

%% 2-3 step: estimate InverseGaussian distribution of stim_on epochs and exclude bad epochs where regressed power exceed the distribution
%-- Compute max power over time for each epoch
max_pows_conf = max_pows_stim_on;       max_pows_conf(outlier_idx) = nan;

%-- Fit distribution
%%% fit to InverseGaussian distribution
for el = 1:size(max_pows_conf,2)
    estdist = fitdist(max_pows_conf(:,el),'InverseGaussian');
    outlier_thresh(el) = estdist.icdf(1-threshout);
end
% %%% fit to empirical distribution
% for el = 1:size(max_pows_conf,2)
%     [estdist,estpoint] = ecdf(max_pows_conf(:,el));
%     outlier_thresh(el) = estpoint(find(estdist>(1-threshdist),1));
% end

outlier_idx2 = max_pows > outlier_thresh;

%% Output
outlier_idx = or(outlier_idx1,outlier_idx2);
epochs(:,outlier_idx) = nan;


end
