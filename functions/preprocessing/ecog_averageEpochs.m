function [epochs_averaged, epochs_sd, stim_idx] = ecog_averageEpochs(epochs, events, stimIDs, avgrange)   

% Computes an average and standard deviations of epoched ecog timecourses
% for a set of for stimulus conditions by matching stimIDs to a BIDS-style
% events table
%
% [epochs_averaged, epochs_sd, stim_idx] = ecog_averageEpochs(epochs, events, stimIDs, [avgrange])   
% 
% INPUT (required)
% epochs  : matrix array of dimensions time x events x channels
% events  : BIDS-formatted events table
% stimIDs : cell array containing either a list of numbers indicating
%           events.trial_types or names indicating events.trial_names
% 
% INPUT (option)
% avgrange : name of range in which data is averaged
%            'all'(default), 'sessions', 'runs'
%            events.session_name and events.run_name are required
%            output data is arranged as stimuli x groups
%
% see also ecog_makeEpochs.m and ecog_normalizeEpochs.m

if ~iscell(stimIDs), stimIDs = {stimIDs}; end
if ~exist('avgrange','var') || isempty(avgrange), avgrange = 'all'; end

switch lower(avgrange)
    case {'all'},       rpttbl = ones(height(events),1);
    case {'sessions'},  rpttbl = findgroups(events.session_name);
    case {'runs'},      rpttbl = findgroups(events(:,{'session_name','run_name'}));
    otherwise,          error('''%s'' is invalid.',avgrange);
end

[nsamp, ~, nchan] = size(epochs);
nstim = length(stimIDs);
nrpt  = length(unique(rpttbl));
navg  = nstim * nrpt;
epochs_averaged = nan(nsamp, navg, nchan);
epochs_sd       = nan(nsamp, navg, nchan);
stim_idx = cell(navg,1);

for ii = 1:navg
    istim = mod(ii-1,nstim)+1;
    irpt  = ceil(ii./nstim);
    if ~isnumeric(stimIDs)
        trial_idx = contains(events.trial_name, stimIDs{istim});        
    else
        trial_idx = events.trial_type == stim_names(istim);
    end
    trial_idx = trial_idx & rpttbl==irpt;
    epochs_averaged(:,ii,:) = mean(epochs(:,trial_idx,:),2, 'omitnan');
    epochs_sd(:,ii,:) = std(epochs(:,trial_idx,:),0,2,'omitnan');
    stim_idx{ii} = find(trial_idx);
end

end

