function [epochs, channels, chan_idx, R2, epochs_split] = ecog_selectElectrodes(epochs, channels, events, t, opts)
% Electrode selection
%
% Select electrodes based on either a simple threshold (cf Zhou et al,
% 2019), a split half method, or single-trial variance relative to mean (cf
% Stigliani et al., 2019). 
% Note: each option requires different combinations of inputs and fields of opts to be set.
% - thresh: requires t, opts.elec_max_thresh, opts.elec_mean_thresh
% - splithalf: requires events, opts.stimnames, opts.elec_splithalf_thresh
% - meanpredict: requires events, opts.stimnames, opts.elec_meanpredict_thresh

if ~exist('events','var'), events = []; end
if ~exist('t','var'), t = []; end
if ~exist('opts','var'), opts = []; end
    
method = opts.elec_selection_method;

switch method
    
    case 'thresh'
        
        % Requires t, opts.elec_max_thresh and opts.elec_mean_thresh to be defined
        
        % Check:
        if ~exist('t', 'var') || isempty(t)
            error('Electrode selection method %s requires t as input', method);
        end
        if ~isfield(opts, 'elec_max_thresh') || ~isfield(opts, 'elec_mean_thresh')
            error('Electrode selection method %s requires definition of max and mean thresholds', method);
        end
        
        % Compute mean across all trials
        mean_resp = mean(epochs,2,'omitnan');

        % Initialize selection to include all channels
        chan_idx = ones(height(channels),2);

        % Exclude channels based on max and mean:
        stim_on_idx = t > opts.stim_on(1) & t <= opts.stim_on(2);  
        chan_idx(:,1) = max(mean_resp(stim_on_idx,:),[],1) > opts.elec_max_thresh;
        chan_idx(:,2) = mean(mean_resp(stim_on_idx,:),1) > opts.elec_mean_thresh;

        % Combine criteria
        chan_idx = sum(chan_idx,2) == size(chan_idx,2);
        
        % Add nans to channel table for consistency with other methods
        channels.noiseceilingR2 = nan(height(channels),1);
        epochs_split = [];
    
    case 'splithalf'
        
        % Requires events, opts.stimnames, opts.elec_splithalf_thresh
        
        % Check:
        if ~exist('events', 'var') || isempty(events)
            error('Electrode selection method %s requires events table', method);
        end
        if ~isfield(opts, 'stimnames') || isempty(opts.stimnames)
            error('Electrode selection method %s requires definition of stimnames', method);
        end
        if ~isfield(opts, 'elec_splithalf_thresh') 
            error('Electrode selection method %s requires definition of splithalf threshold', method);
        end     
        
        % Get trial indices for each stimname
        [~,~,trial_idx] = ecog_averageEpochs(epochs, events, opts.stimnames);
        [nSamp, ~, nChan] = size(epochs);
        
        % Put trials in the last dimension
        temp_epochs = permute(epochs,[3 1 2]);
        
        epochs_split = nan([2 size(temp_epochs,1) size(temp_epochs,2) length(opts.stimnames)]);
        % Average first and second halfs of trials for each stimulusname
        for stim = 1:length(opts.stimnames)
            idx = trial_idx{stim};
            % This way of splitting the data will mix runs and sessions. I
            % think that's OK given that we convert trials to percent
            % signal change based on baselines computed within runs first.
            trial_idx1 = idx(1:2:length(idx));
            trial_idx2 = idx(2:2:length(idx));
            epochs_split(1,:,:,stim) = mean(temp_epochs(:,:,trial_idx1),3,'omitnan');
            epochs_split(2,:,:,stim) = mean(temp_epochs(:,:,trial_idx2),3,'omitnan');
        end
        
        % Correlate across halfs

        % X1: chans x concatenated stim half 1
        X1 = squeeze(epochs_split(1,:,:));
        % X2: chans x concatenated stim half 2
        X2 = squeeze(epochs_split(2,:,:));
        
        % % Pairwise correlation between all channels across the two sets
        %[r] = corr(X1',X2');
        % % Select diagonal (correlation of channel with itself)
        %R2 = r(eye(size(r))==1);
        
        % Compute R2
        R2x = computeR2(X1',X2');
        R2y = computeR2(X2',X1');
        R2 = mean([R2x' R2y'], 2);
        % KK code:
        %R2x = calccod(X1', X2');
        %R2y = calccod(X1', X2');
        %R2 = mean([R2x' R2y'], 2);

        % Add values to channel table
        channels.noiseceilingR2 = round(R2,2);
        
        % Select channels based on threshold
        chan_idx = R2 >= opts.elec_splithalf_thresh;
        
%         % Debug       
%         for el = 1:nChan
%             figure;
%             figureName = 'compareR2andPearson';
%             subplot(2,1,1); bar(r); title('Pearson correlation'); xlabel('channel');
%             set(gca, 'XTick', 1:nChan, 'XTickLabel', channels.name);
%             subplot(2,1,2); bar(R2); title('R2'); xlabel('channel');
%             set(gca, 'XTick', 1:nChan, 'XTickLabel', channels.name);
%             set(gcf, 'Position', get(0, 'Screensize'));
%             set(findall(gcf,'-property','FontSize'),'FontSize',14)
%             saveas(gcf, fullfile(plotSaveDir, figureName), 'png'); close;
%         end
        
	case 'meanpredict'
        
        % Requires events, opts.stimnames, opts.elec_meanpredict_thresh
        
        % Check:
        if ~exist('events', 'var') || isempty(events)
            error('Electrode selection method %s requires events table', method);
        end
        if ~isfield(opts, 'stimnames') || isempty(opts.stimnames)
            error('Electrode selection method %s requires definition of stimnames', method);
        end
        if ~isfield(opts, 'elec_meanpredict_thresh') 
            error('Electrode selection method %s requires definition of meanpredict threshold', method);
        end     
        
        % Get trial indices for each stimname
        [~,~,trial_idx] = ecog_averageEpochs(epochs, events, opts.stimnames);
        [~, ~, nChan] = size(epochs);
        nStim = length(opts.stimnames);
        R2 = nan(nChan,1);
        % Compute 'noise ceiling' for each electrode as per the method
        % described in Stigliani et al., 2017, PNAS (code adapted from
        % https://github.com/VPNL/TemporalChannels, tchROI.m)
        
        for ee = 1:nChan
            
            % Find average response for each trial type
            % Use cell arrays in case there aren't the same number of trials
            % for each condition            
            trials = cell(nStim,1); nTrials = cell(nStim,1);
            for stim = 1:length(opts.stimnames)
                idx = trial_idx{stim};
                trials{stim,1} = epochs(:,idx,ee);
                nTrials{stim} = length(idx);
            end
            % Find average response for each trial type
            trials_avg = cellfun(@(X) mean(X, 2, 'omitnan'), trials, 'uni', false);
            trials_avg_rep = cellfun(@(X,y) repmat(X,[1 y]), trials_avg, nTrials, 'uni', false);
            % Calculate residual between individual trials and means
            trials_err = cellfun(@(X, Y) X - Y, trials, trials_avg_rep, 'uni', false);
            trials_err = cellfun(@(X) sum(sum(X .^ 2, 'omitnan')), trials_err, 'uni', false);
            % Calculate variance explained by mean trial responses
            total_err = sum(cell2mat(trials_err));
            trial_mean = mean(mean(cell2mat(trials), 'omitnan'));
            trials_var = cell2mat(trials) - trial_mean;
            total_var = sum(sum(trials_var .^ 2,'omitnan'));
            R2(ee) = 1 - (total_err / total_var);

        end
        
        % Select channels based on threshold
        chan_idx = R2 >= opts.elec_meanpredict_thresh;
        
        % Add values to channel table
        channels.noiseceilingR2 = round(R2,2);
        epochs_split = [];
end

% Exclude depth electrodes
%if opts.elec_exclude_depth 
%    chan_idx(contains(lower(channels.type), 'seeg')) = 0;
%end  

% Apply chan_idx to epochs and update channels table
epochs = epochs(:,:,chan_idx);
channels = channels(chan_idx,:);

end