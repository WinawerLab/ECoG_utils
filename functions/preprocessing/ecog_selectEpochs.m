function [select_idx, max_epochs, thresh] = ecog_selectEpochs(epochs_v, epochs_b, t, opts)
% Epoch selection 
%    
% Remove epochs whose max response amplitude is opts.epoch_outlier_thresh x more
% or less than the median response in that channel

% Two exclusion criteria:

% 1. Remove epochs in which there is a sudden jump in voltage from one
% timepoint to the next.

if ~isfield(opts, 'epoch_jump_thresh')  || isempty(opts.epoch_jump_thresh)
    % do not apply this criterion
    warning('[%s] No jump threshold defined, not detecting jumps.', mfilename);

    outlier_idx = zeros(size(epochs_v));

else
    [nSamp, nStim, nChans] = size(epochs_v);

    temp_epochs = reshape(epochs_v, [nSamp nStim*nChans]);
    diff_epochs = abs(diff(temp_epochs, 1, 1));
    outlier_idx = any(diff_epochs > opts.epoch_jump_thresh, 1);
    outlier_idx = reshape(outlier_idx, [nStim nChans]);         
end

epochs_b(:,outlier_idx) = nan;

% 2. Remove epochs in which the absolute maximum amplitude is in the tail
% of the distribution 

if ~isfield(opts, 'epoch_outlier_thresh')  || isempty(opts.epoch_outlier_thresh)
    % do not apply this criterion
    warning('[%s] No outlier threshold defined, not detecting outlier epochs.', mfilename);
    thresh = []; max_epochs = [];
else

    if ~isfield(opts, 'stim_on') 
        warning('[%s] Epoch selection based on threshold requires definition of stim_on period', filename);
        return
    else
        % Compute time indices
        stim_on_idx = t > opts.stim_on(1) & t <= opts.stim_on(2);  

        % Compute max over time for each epoch
        max_epochs = squeeze(max(epochs_b(~stim_on_idx,:,:),[],1, 'omitnan'));

        % Compute mean over stim_on period across epochs
        max_epochs_stim_on = squeeze(max(epochs_b(stim_on_idx,:,:),[],1, 'omitnan'));
        
        % Put all epochs with max higher than outlier_thresh * mean to nans
        thresh = opts.epoch_outlier_thresh * std(max_epochs_stim_on,[],1, 'omitnan') + mean(max_epochs_stim_on, 'omitnan');
        outlier_idx2 = max_epochs > thresh | max_epochs_stim_on > thresh;             
        %outlier_idx2 = max_epochs > thresh; 
        outlier_idx(outlier_idx2) = 1;
    end
end

% Print number of rejected epochs
nReject = length(find(outlier_idx));
nTotal =  numel(outlier_idx);
percReject = (nReject/nTotal)*100;
fprintf('[%s] Number of removed epochs %d out of %d = %0.3f percent rejections \n',mfilename, nReject, nTotal, percReject);

select_idx = ~outlier_idx;

% DEBUG
% epochs_b2 = epochs_b;
% epochs_b2(:,outlier_idx) = nan;
% for ii = 1:size(epochs_b,3)
%     figure;hold on
%     subplot(1,2,1);
%     plot(epochs_b(:,:,ii));
%     title('original');
%     subplot(1,2,2)
%     plot(epochs_b2(:,:,ii));
%     title('selected');
% end

end


