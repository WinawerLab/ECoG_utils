function [epochs, channels, chan_idx] = ecog_selectElectrodes(epochs, channels, events, t, opts, plotSaveDir)
    
switch opts.elec_selection_method
    
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
        
    case 'splithalf'
        
        % Requires opts.elec_split_thresh 
        % Check:
        if ~exist('events', 'var') || isempty(events)
            error('Electrode selection method %s requires events table', method);
        end
        if ~isfield(opts, 'stimnames') || isempty(opts.stimnames)
            error('Electrode selection method %s requires definition of stimnames', method);
        end
        if ~isfield(opts, 'elec_split_thresh') 
            error('Electrode selection method %s requires definition of split threshold', method);
        end     
        
        % Get trial indices for each stimname
        [~,~,trial_idx] = ecog_averageEpochs(epochs, events, opts.stimnames);
        [nSamp, ~, nChan] = size(epochs);
        
        temp_epochs = permute(epochs,[3 1 2]);
        
        % Average first and second halfs of trials for each stimulusname
        for stim = 1:length(opts.stimnames)
            idx = trial_idx{stim};
            trial_idx1 = idx(1:2:length(idx));
            trial_idx2 = idx(2:2:length(idx));
            epochs_averaged(1,:,:,stim) = mean(temp_epochs(:,:,trial_idx1),3,'omitnan');
            epochs_averaged(2,:,:,stim) = mean(temp_epochs(:,:,trial_idx2),3,'omitnan');
        end
        
        % Correlate across halfs

        % X1: chans x concatenated stim half 1
        X1 = squeeze(epochs_averaged(1,:,:));
        % X2: chans x concatenated stim half 2
        X2 = squeeze(epochs_averaged(2,:,:));
        
        % Pairwise correlation between all channels across the two sets
        [r] = corr(X1',X2');
        
        % Select diagonal (correlation of channel with itself)
        r = r(eye(size(r))==1);
        
        % Compute R2
        R2 = computeR2(X1',X2');
        
        chan_idx = R2 > opts.elec_split_thresh;
        
        for el = 1:nChan
            figureName = sprintf('%s_%s_%s', ...
                channels.name{el}, channels.bensonarea{el}, channels.wangarea{el});
            figure;hold on;
            plot(X1(el,:), 'r','LineWidth', 2);
            plot(X2(el,:), 'b','LineWidth', 2);
            axis tight
            set(gca, 'XTick', 1:nSamp:size(X1,2), 'XTickLabel', opts.stimnames);
            xtickangle(45)
            title(sprintf('%s %s %s Pearson = %0.2f R2 = %0.2f', ...
                channels.name{el}, channels.bensonarea{el}, channels.wangarea{el}, ...
                r(el), R2(el)));
            set(gcf, 'Position', get(0, 'Screensize'));
            set(findall(gcf,'-property','FontSize'),'FontSize',14)
            saveas(gcf, fullfile(plotSaveDir, figureName), 'png'); close;
        end
        
        figure;
        figureName = 'compareR2andPearson';
        subplot(2,1,1); bar(r); title('Pearson correlation'); xlabel('channel');
        set(gca, 'XTick', 1:nChan, 'XTickLabel', channels.name);
        subplot(2,1,2); bar(R2); title('R2'); xlabel('channel');
        set(gca, 'XTick', 1:nChan, 'XTickLabel', channels.name);
        set(gcf, 'Position', get(0, 'Screensize'));
        set(findall(gcf,'-property','FontSize'),'FontSize',14)
        saveas(gcf, fullfile(plotSaveDir, figureName), 'png'); close;
        
	case 'meanpredict'
        
end

% Exclude depth electrodes
if opts.elec_exclude_depth 
    chan_idx(contains(lower(channels.type), 'seeg')) = 0;
end  

% Apply chan_idx to epochs and update channels table
epochs = epochs(:,:,chan_idx);
channels = channels(chan_idx,:);

end