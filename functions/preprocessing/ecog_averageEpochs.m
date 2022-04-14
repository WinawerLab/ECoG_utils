function [epochs_averaged, epochs_sd, stim_idx] = ecog_averageEpochs(epochs, events, stimIDs)   

% Computes an average and standard deviations of epoched ecog timecourses
% for a set of for stimulus conditions by matching stimIDs to a BIDS-style
% events table
%
% [epochs_averaged, epochs_sd, stim_idx] = ecog_averageEpochs(epochs, events, stimIDs)   
% 
% INPUT (required)
% epochs  : matrix array of dimensions time x events x channels
% events  : BIDS-formatted events table
% stimIDs : cell array containing either a list of numbers indicating
%           events.trial_types or names indicating events.trial_names
%
% see also ecog_makeEpochs.m and ecog_normalizeEpochs.m

if ~iscell(stimIDs), stimIDs = {stimIDs}; end

[nsamp, ~, nchan] = size(epochs);
nstim = length(stimIDs);
epochs_averaged = nan(nsamp, nstim, nchan);

epochs_sd = nan(size(epochs,1), length(stimIDs), size(epochs,3));
stim_idx = cell(length(stimIDs),1);

for ii = 1:length(stimIDs)
    if ~isnumeric(stimIDs)
        trial_idx = find(contains(events.trial_name, stimIDs{ii}));        
    else
        trial_idx = find(events.trial_type == stim_names(ii));
    end
    epochs_averaged(:,ii,:) = mean(epochs(:,trial_idx,:),2, 'omitnan');
    epochs_sd(:,ii,:) = std(epochs(:,trial_idx,:),0,2,'omitnan');
    stim_idx{ii} = trial_idx;
end

end

