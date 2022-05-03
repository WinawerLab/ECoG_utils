function [chan_idx, R2, epochs_split] = ecog_selectElectrodes(epochs, channels, events, t, opts)
% Electrode selection
%
% Select electrodes based on either a simple threshold (cf Zhou et al,
% 2019), a split half method, or single-trial variance relative to mean (cf
% Stigliani et al., 2019). 
% Note: each option requires different combinations of inputs and fields of opts to be set.
% - thresh: requires t, opts.stim_on, opts.elec_max_thresh, opts.elec_mean_thresh
% - splithalf: requires events, opts.stimnames, opts.elec_splithalf_thresh, opts.elec_splitmethod
% - meanpredict: requires events, opts.stimnames, opts.elec_meanpredict_thresh

narginchk(4,5);
if ~exist('opts','var'), opts = []; end
if nargin == 4 && isstruct(t)
    opts    = t;
    t       = [];
end
assert(~isempty(opts) && isfield(opts,'elec_selection_method') && ~isempty(opts.elec_selection_method),...
    'opts.elec_selection_method is required.');
method = opts.elec_selection_method;
R2 = [];
epochs_split = [];

switch method
    
    case 'thresh'
        
        % Requires t, opts.elec_max_thresh and opts.elec_mean_thresh to be defined
        
        % Check:
        if isempty(t)
            if ~isepmty(events) && ~istable(events)
                t       = events;
            else
                error('Electrode selection method %s requires t as input', method);
            end
        end
        if ~isfield(opts, 'stim_on') || isempty(opts.stim_on)
            Warning('Electrode selection method %s requires definition of stim_on', method);
            opts.stim_on = [min(t)-1 max(t)];
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
    
    case 'splithalf'
        
        % Requires events, opts.stimnames, opts.elec_splithalf_thresh, opts.elec_splitmethod
        
        % Check:
        if isempty(events)
            error('Electrode selection method %s requires events table', method);
        end
        if ~isfield(opts, 'stimnames') || isempty(opts.stimnames)
            warning('Electrode selection method %s requires definition of stimnames', method);
            opts.stimnames = unique(events.trial_name);
        end
        if ~isfield(opts, 'elec_splithalf_thresh') 
            error('Electrode selection method %s requires definition of splithalf threshold', method);
        end     
        if ~isfield(opts, 'elec_splitmethod') 
            opts.elec_splitmethod = 'alternative';
        end     
        
        % Get trial indices for each stimname
        elec_splitmethod = lower(opts.elec_splitmethod);
        if endsWith(elec_splitmethod,{'inruns','insessions'})
            elec_splitmethod = strsplit(elec_splitmethod,'in');
            avgrange = elec_splitmethod{end};
            elec_splitmethod = elec_splitmethod{1};
        else
            avgrange = 'all';
        end
        [~,~,trial_idx] = ecog_averageEpochs(epochs, events, opts.stimnames, avgrange);
        nstims = length(trial_idx);
        
        % Put trials in the last dimension
        temp_epochs = permute(epochs,[3 1 2]);
        
        epochs_split = nan([2 size(temp_epochs,1) size(temp_epochs,2) nstims]);
        % Average first and second halfs of trials for each stimulusname
        for stim = 1:nstims
            idx = trial_idx{stim};
            % This way of splitting the data will mix runs and sessions. I
            % think that's OK given that we convert trials to percent
            % signal change based on baselines computed within runs first.
            switch elec_splitmethod
                case {'alternative'}
                    trial_idx1 = idx(1:2:length(idx));
                    trial_idx2 = idx(2:2:length(idx));
                case {'twohalves'}
                    trial_idx1 = idx(1:fix(length(idx)./2));
                    trial_idx2 = idx(ceil(length(idx)./2):length(idx));
                otherwise
                    error('%s is unknown split method',opts.elec_splitmethod);
            end
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
        if isempty(events)
            error('Electrode selection method %s requires events table', method);
        end
        if ~isfield(opts, 'stimnames') || isempty(opts.stimnames)
            warning('Electrode selection method %s requires definition of stimnames', method);
            opts.stimnames = unique(events.trial_name);
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
            trials_err = cellfun(@(X) sum(sum(X .^ 2, 'omitnan'), 'omitnan'), trials_err, 'uni', false);
            % Calculate variance explained by mean trial responses
            total_err = sum(cell2mat(trials_err));
            trial_mean = mean(mean(cell2mat(trials'), 'omitnan'), 'omitnan');
            trials_var = cell2mat(trials') - trial_mean;
            total_var = sum(sum(trials_var .^ 2,'omitnan'), 'omitnan');
            R2(ee) = 1 - (total_err / total_var);

        end
        
        % Select channels based on threshold
        chan_idx = R2 >= opts.elec_meanpredict_thresh;
        
    otherwise
        error('%s is unknown method.',method);

end

% Print number of rejected electrodes
nReject = length(find(~chan_idx));
nTotal =  numel(chan_idx);
percReject = (nReject/nTotal)*100;
fprintf('[%s] Number of removed electrodes %d out of %d = %0.1f percent rejections \n',mfilename, nReject, nTotal, percReject);
fprintf('[%s] Number of kept electrodes %d out of %d = %0.1f percent kept \n',mfilename, nTotal-nReject, nTotal, 100-percReject);

end