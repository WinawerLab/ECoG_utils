function [epochs, outlier_idx, max_epochs, thresh] = ecog_selectEpochs(epochs, t, opts)
% Epoch selection 
%    
% Remove epochs whose max response amplitude is thresh x more
% or less than the median response in that channel

if ~isfield(opts, 'stim_on') 
    warning('[%s] Epoch selection requires definition of stim_on period, not performing epoch selection', filename);
    return
    
else

    % Compute time indices
    stim_on_idx = t > opts.stim_on(1) & t <= opts.stim_on(2);  

    % Two exclusion criteria:

    % 1. Remove epochs in which there is a sudden jump in voltage from one
    % timepoint to the next.
    if ~isfield(opts, 'epoch_jump_thresh')  || isempty(opts.epoch_jump_thresh)
        % do not apply this criterion
        warning('[%s] No jump threshold defined, not detecting jumps.', mfilename);

        outlier_idx = ones(size(epochs));
    
    else
        [nSamp, nStim, nChans] = size(epochs);

        temp_epochs = reshape(epochs, [nSamp nStim*nChans]);
        diff_epochs = abs(diff(temp_epochs, 1, 1));
        outlier_idx = any(diff_epochs > opts.epoch_jump_thresh, 1);
        outlier_idx = reshape(outlier_idx, [nStim nChans]);         
    end

    % 2. Remove epochs in which the absolute maximum amplitude is in the tail
    % of the distribution (cg KY)

    if ~isfield(opts, 'epoch_outlier_thresh')  || isempty(opts.epoch_outlier_thresh)
        % do not apply this criterion
        warning('[%s] No outlier threshold defined, not detecting outlier epochs.', mfilename);

    else
        
        % Compute max over stim_on period across epochs
        max_epochs = squeeze(max(epochs,[],1));

        % Compute max over stim_on period across epochs
        max_epochs_stim_on = squeeze(max(epochs(stim_on_idx,:,:),[],1));
        mn_max = mean(max_epochs_stim_on,1);
        % Compute difference between epochs max and the average max
        % amplitude in the stim on period
        
        % Put all epochs with max higher than outlier_thresh * median to nans
        thresh = prctile(abs(max_epochs(:)), opts.epoch_outlier_thresh);
        outlier_idx2 = abs(max_epochs) > thresh;             
        outlier_idx(outlier_idx2) = 1;
    end

    epochs(:,outlier_idx) = nan;
end

end
