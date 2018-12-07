function [out] = ecog_plotTrials(trials, whichElectrodes, trialType, collapseTrialTypes, smoothingLevelInMs, baselineType)

if nargin < 6 || isempty(baselineType)
    baselineType = 'all';
end
if nargin < 5 || isempty(smoothingLevelInMs)
    smoothLevel = 0;
else 
    smoothLevel = round(smoothingLevelInMs/(1000/trials.fsample)); 
    % number of samples to smooth over
end

if nargin < 4 || isempty(collapseTrialTypes)
    collapseTrialTypes = 'no';
end

% Find electrodes in data
el_index = ecog_matchChannels(whichElectrodes, trials);

trial_index = [];
% Find trials in data
if isempty(trialType)
    trialType = {'all trials'};
	trial_index{1} = 1:size(trials.events,1);
    %baseline_index = vertcat(trial_index{:});
    baseline_index = find(~contains(trials.events.trial_name, 'PRF'));
else
    for ii = 1:length(trialType)
        trial_index{ii} = find(contains(trials.events.trial_name, trialType{ii}));
    end
    %fprintf('[%s] Matching whichTrials to events.trial_name\n', mfilename)    
    baseline_index = vertcat(trial_index{:});
        
    switch collapseTrialTypes
        case 'yes'
            trial_index = [];
            trial_index{1} = baseline_index;
    end 
end

% Plot specs
colors = copper(length(trial_index));

out = struct;

% Plot
for dataType = {'broadband', 'evoked'}
    
    thisDataType = dataType{:};
    
    figure('Name', thisDataType); 
    
    for ii = 1:length(el_index)

        % Make a separate plot for each channel
        subplot(ceil(sqrt(length(el_index))),ceil(sqrt(length(el_index))), ii); hold on;
        %subplot(ceil(sqrt(length(el_index))),floor(sqrt(length(el_index))), ii); hold on;
        
        elData = squeeze(trials.(thisDataType)(el_index(ii),:,:));
        
        % Baseline correction: 
        switch baselineType
            case 'all'
                baseline_index = find(~contains(trials.events.trial_name, 'PRF'));
                baseline = mean(mean(elData(trials.time<0,baseline_index),1),2);
            case 'selectedtrials'
                baseline = mean(mean(elData(trials.time<0,baseline_index),1),2);
        end
        
        switch thisDataType
            case 'broadband'
                % percent signal? (similar to 'relchange' in fieldtrip)
                elData = (elData - baseline) ./ baseline;
            case 'evoked'
                % standard ERP approach?
                elData = (elData - baseline);
        end

        % Select subset of trials to plot
        for jj = 1:length(trial_index)
            % Compute mean and standard error of the mean
            mnToPlot(:,jj) = median(elData(:,trial_index{jj}),2);
            Llim(:,jj) = mnToPlot(:,jj)-(std(elData(:,trial_index{jj}),0,2)/sqrt(size(elData(:,trial_index{jj}),2)));
            Ulim(:,jj) = mnToPlot(:,jj)+(std(elData(:,trial_index{jj}),0,2)/sqrt(size(elData(:,trial_index{jj}),2)));    
            
        end

        % Smooth the data?
        if smoothLevel > 0
            temp = smooth(mnToPlot, smoothLevel); mnToPlot = reshape(temp, size(mnToPlot));
            temp = smooth(Llim, smoothLevel); Llim = reshape(temp, size(Llim));
            temp = smooth(Ulim, smoothLevel); Ulim = reshape(temp, size(Ulim));
        end
        
        % Plot standard errors
        for jj = 1:length(trial_index)
            h = ciplot(Llim(:,jj),Ulim(:,jj),trials.time,colors(jj,:), 0.25);
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        end
        
        % Plot means
        for jj = 1:length(trial_index)
            plot(trials.time, mnToPlot(:,jj),'Color', colors(jj,:), 'LineWidth',2);
        end
        if ii == 1
            xlabel('time (s)');
            switch thisDataType
                case 'broadband'
                    ylabel(['relative change in broadband power ( ' num2str(min(trials.bb_bands(:))) '-' num2str(max(trials.bb_bands(:))) ' Hz)']);
                case 'evoked'
                    ylabel('evoked response amplitude');
            end
        end

        % Check if these electrodes have matches with visual atlases, if so, add
        % that to the plot title
         if iscell(whichElectrodes)
            electrodeName = whichElectrodes{ii};
        else
            electrodeName = whichElectrodes;
        end
        viselec_name = [];
        for atlas = {'wang2015_atlas','benson14_varea', 'wang15_mplbl'}
            if ~isempty(trials.viselec)
                if ~isempty(trials.viselec.(atlas{:}))
                    %viselec = contains(trials.viselec.(atlas{:}).elec_labels, electrodeName);
                    viselec = find(strcmp(electrodeName,trials.viselec.(atlas{:}).elec_labels));
                    if any(viselec)
                        %viselec_name = [viselec_name atlas{:}(1:8) ':' trials.viselec.(atlas{:}).area_labels{viselec} ' '];
                        atlasstr = strsplit(atlas{:},'_');
                        viselec_name = [viselec_name atlasstr{1} ':' trials.viselec.(atlas{:}).area_labels{viselec} ' '];
                    end
                end
            end
        end
       
        title([electrodeName ' ' viselec_name]);

        % Set y-axis limits
        lim = [min(mnToPlot(:)) max(mnToPlot(:))];
        ylim = [lim(1)-(0.2*lim(1)*sign(lim(1))) lim(2)+(0.2*lim(2)*sign(lim(2)))];
        set(gca, 'YLim', ylim);

        % Add stim onset and zero lines
        line([0 0], ylim,'LineStyle', ':', 'Color', 'k');
        line([trials.time(1) trials.time(end)], [0 0],'LineStyle', ':', 'Color', 'k');

        % Add legend
        if ii == 1
            switch collapseTrialTypes
                case 'no'
                    legend(trialType);
                case 'yes'
                    legend([trialType{:}])
            end
        end
        
        out.(thisDataType).(whichElectrodes{ii}).mn = mnToPlot';
        out.(thisDataType).(whichElectrodes{ii}).se = (Ulim-mnToPlot)';
        set(gca, 'XLim', [-0.5 1.5]);
        set(gca, 'FontSize', 12);

    end
    %set(gcf, 'Position', [150 100 1500 1250]);
    set(gcf, 'Position', [150 100 750 625]);
	out.time = trials.time;
    
end

% For broadband, add figure with sum across time:

% define time window to sum over: 0-1sec?

% figure('Name', thisDataType); 
%     for ii = 1:length(el_index)
%     end
    
