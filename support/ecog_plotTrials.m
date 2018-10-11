function [out] = ecog_plotTrials(trials, whichElectrodes, trialType, collapseTrialTypes, smoothLevel)

if nargin < 5 || isempty(smoothLevel)
    smoothLevel = 0;
else 
    smoothLevel = ceil(smoothLevel);
end

if nargin < 4 || isempty(collapseTrialTypes)
    collapseTrialTypes = 'no';
end

% Find electrodes in data
el_index = ecog_matchChannels(whichElectrodes, trials);

trialindex = [];
% Find trials in data
if isempty(trialType)
	trial_index = 1:size(trials.events,1);
else
    trial_index = [];
    switch collapseTrialTypes
        case 'yes'
            trial_index{1} = (contains(trials.events.trial_name, trialType));
        case 'no'
            for ii = 1:length(trialType)
                trial_index{ii} = find(contains(trials.events.trial_name, trialType{ii}));
            end
    end 
end

% Plot specs
colors = copper(length(trial_index));

out = struct;

% Plot
for dataType = {'broadband'}%, 'evoked'}
    
    thisDataType = dataType{:};
    figure('Name', thisDataType); 
    
    for ii = 1:length(el_index)

        % Make a separate plot for each channel
        subplot(ceil(sqrt(length(el_index))),ceil(sqrt(length(el_index))), ii); hold on;

        % Baseline correction: across all trials 
        elData = squeeze(trials.(thisDataType)(el_index(ii),:,:));
        baseline = mean(mean(elData(trials.time<0,:),1),2);
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
            mnToPlot(:,jj) = mean(elData(:,trial_index{jj}),2);
            llim(:,jj) = mean(elData(:,trial_index{jj}),2)-std(elData(:,trial_index{jj}),0,2)/sqrt(size(elData(:,trial_index{jj}),2));
            ulim(:,jj) = mean(elData(:,trial_index{jj}),2)+std(elData(:,trial_index{jj}),0,2)/sqrt(size(elData(:,trial_index{jj}),2));    
        end

        % Smooth the data?
        if smoothLevel > 0
            temp = smooth(mnToPlot, smoothLevel); mnToPlot = reshape(temp, size(mnToPlot));
            temp = smooth(llim, smoothLevel); llim = reshape(temp, size(llim));
            temp = smooth(ulim, smoothLevel); ulim = reshape(temp, size(ulim));
        end
        
        % Plot standard errors
        for jj = 1:length(trial_index)
            h = ciplot(llim(:,jj),ulim(:,jj),trials.time,colors(jj,:), 0.25);
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        end
        
        % Plot means
        for jj = 1:length(trial_index)
            plot(trials.time, mnToPlot(:,jj),'Color', colors(jj,:), 'LineWidth',2);
        end
        xlabel('time (s)');
        switch thisDataType
            case 'broadband'
                ylabel(['relative change in broadband power ( ' num2str(min(trials.bb_bands(:))) '-' num2str(max(trials.bb_bands(:))) ' Hz)']);
            case 'evoked'
                ylabel('evoked response amplitude');
        end
        
        % Check if these electrodes have matches with visual atlases, if so, add
        % that to the plot title
        viselec_name = [];
        for atlas = {'wang2015_atlas','benson14_varea'}
            viselec = contains(trials.viselec.(atlas{:}).elec_labels, whichElectrodes(ii));
            if any(viselec)
                viselec_name = [viselec_name atlas{:}(1:8) ':' trials.viselec.(atlas{:}).area_labels{viselec} ' '];
            end
        end
        title([whichElectrodes{ii} ' ' viselec_name]);

        % Set y-axis limits
        lim = [min(mnToPlot(:)) max(mnToPlot(:))];
        ylim = [lim(1)-(0.2*lim(1)*sign(lim(1))) lim(2)+(0.2*lim(2)*sign(lim(2)))];
        set(gca, 'YLim', ylim);

        % Add stim onset and zero lines
        line([0 0], ylim,'LineStyle', ':', 'Color', 'k');
        line([trials.time(1) trials.time(end)], [0 0],'LineStyle', ':', 'Color', 'k');

        % Add legend
        switch collapseTrialTypes
            case 'no'
                legend(trialType);
            case 'yes'
                legend([trialType{:}])
        end
    end
    set(gcf, 'Position', [150 100 1500 1250]);
    
    out.(thisDataType).mn = mnToPlot';
    out.(thisDataType).se = (ulim-mnToPlot)';
   
end

