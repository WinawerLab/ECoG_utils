
function [out] = ecog_plotTimecourses(trials, whichElectrodes, trialType, specs)

%function [out] = ecog_plotTimecourses(trials, whichElectrodes, trialType, collapseTrialTypes, smoothingLevelInMs, baselineType)
% description

if nargin < 4
    specs = [];
end

if ~isfield(specs, 'dataTypes') || isempty(specs.dataTypes)
    specs.dataTypes = {'broadband', 'evoked'};
end

if ~isfield(specs, 'baselineType') || isempty(specs.baselineType)
    specs.baselineType = 'all';
end

if ~isfield(specs, 'smoothingLevelInMs') || isempty(specs.smoothingLevelInMs)
    specs.smoothLevel = 0;
else 
    specs.smoothLevel = round(specs.smoothingLevelInMs/(1000/trials.fsample)); 
    % number of samples to smooth over
end

if ~isfield(specs, 'collapseTrialTypes') || isempty(specs.collapseTrialTypes)
    specs.collapseTrialTypes = 'no';
end

if ~isfield(specs, 'plot') || isempty(specs.plot)
    specs.plot = [];
end

if ~isfield(specs.plot, 'colorMap') || isempty(specs.plot.colorMap), specs.plot.colorMap = 'copper'; end
if ~isfield(specs.plot, 'nSubPlots') || isempty(specs.plot.nSubPlots), specs.plot.nSubPlots = []; end
if ~isfield(specs.plot, 'addEccToTitle') || isempty(specs.plot.addEccToTitle), specs.plot.addEccToTitle = 'no'; end
if ~isfield(specs.plot, 'showMax') || isempty(specs.plot.showMax), specs.plot.showMax = 'no'; end
if ~isfield(specs.plot, 'fontSize') || isempty(specs.plot.fontSize), specs.plot.fontSize = 12; end
if ~isfield(specs.plot, 'XLim') || isempty(specs.plot.XLim), specs.plot.XLim = [-0.2 1];end
if ~isfield(specs.plot, 'YLim'), specs.plot.YLim = [];end
if ~isfield(specs.plot, 'XScale') || isempty(specs.plot.XScale), specs.plot.XScale = 'linear';end
if ~isfield(specs.plot, 'YScale') || isempty(specs.plot.YScale), specs.plot.YScale = 'linear';end

out = struct;

% Find electrodes in data
el_index = ecog_matchChannels(whichElectrodes, trials);
if isempty(el_index)
    fprintf('[%s] Did not find any matching channels, not plotting \n', mfilename); 
    return
end
% Find trials in data
trial_index = [];
if isempty(trialType)
    trialType = {'all'};
	trial_index{1} = 1:size(trials.events,1);
    %baseline_index = vertcat(trial_index{:});
    baseline_index = find(~contains(lower(trials.events.task_name), {'prf'}));
else
    for ii = 1:length(trialType)
        trial_index{ii} = find(contains(trials.events.trial_name, trialType{ii}));
        %trial_index{ii} = find(contains(trials.events.task_name, trialType{ii}));
    end
    %HACK
%     trial_index{1} = 100:114;%trial_index{1}(33:64);
%     trial_index{1} = 85:99;%trial_index{1}(33:64);
%     trial_index{1} = 70:84;%trial_index{1}(33:64);
%     trial_index{1} = 55:69;%trial_index{1}(33:64);
%     fprintf('[%s] Matching whichTrials to events.trial_name\n', mfilename)    
    baseline_index = vertcat(trial_index{:});
        
    switch specs.collapseTrialTypes
        case 'yes'
            trial_index = [];
            trial_index{1} = baseline_index;
    end 
end

% Plot settings
cmap = eval(specs.plot.colorMap);
colors = cmap(1:round((length(cmap)/length(trial_index))):length(cmap),:);
%colors = sortrows(colors,'descend');

out.elecs = whichElectrodes;
out.titles = cell(size(out.elecs));

% Decide how many subplots are needed
if ~isempty(specs.plot.nSubPlots)
    nRow = specs.plot.nSubPlots(1);
    nCol = specs.plot.nSubPlots(2);
else
    nPlot = length(el_index);
    nRow = ceil(sqrt(nPlot));
    nCol = ceil(sqrt(nPlot));
    if nPlot <= (nRow*nCol)-nCol
        nRow = nRow-1;
    end
end

% Plot
for d = 1:length(specs.dataTypes)
    
    thisDataType = specs.dataTypes{d};
    figureName = strsplit(trialType{1},'-');
    figure('Name', [figureName{1} ' ' thisDataType]); 
    
    subplotWithLegendIndex = 1;
    for ii = 1:length(el_index)

        % Make a separate plot for each channel
        subplot(nRow,nCol,ii); hold on;
        %subplot(ceil(sqrt(length(el_index))),floor(sqrt(length(el_index))), ii); hold on;
        %subplot(8,8,ii); hold on;
        
        % Add x and y labels
        if ii == 1
            xlabel('time (s)');
            switch thisDataType
                case 'broadband'
                    ylabel(['change in power ( ' num2str(min(trials.bb_bands(:))) '-' num2str(max(trials.bb_bands(:))) ' Hz)']);
                case 'evoked'
                    ylabel('evoked response amplitude');
            end
        end
        
        % Check if we have data to plot for this electrode, if not, go on
        % to the next one
        if el_index(ii)>0
            elData = squeeze(trials.(thisDataType)(el_index(ii),:,:));
        else
            subplotWithLegendIndex = subplotWithLegendIndex+1;
            set(gca, 'FontSize', specs.plot.fontSize);
            continue
        end
        
        % Baseline correction: 
        switch specs.baselineType
            case 'all'
                baseline_index = find(~contains(trials.events.trial_name, {'PRF', 'BLANK'}));
                baseline = nanmean(nanmean(elData(trials.time<0,baseline_index),1),2);
            case 'selectedtrials'
                baseline = nanmean(nanmean(elData(trials.time<0,baseline_index),1),2);
        end
        
        switch thisDataType
            case 'broadband'
                % percent signal? (similar to 'relchange' in fieldtrip)
                %elData = (elData - baseline) ./ baseline;
                elData = (elData - nanmean(elData(trials.time<0,:),1)) ./ baseline;
            case 'evoked'
                % standard ERP approach(?)
                %elData = (elData - baseline);
                elData = elData - nanmean(elData(trials.time<0,:),1);
        end

        % Select subset of trials to plot
        for jj = 1:length(trial_index)
            % Compute mean and standard error of the mean
            mnToPlot(:,jj) = nanmean(elData(:,trial_index{jj}),2);
            Llim(:,jj) = mnToPlot(:,jj)-(nanstd(elData(:,trial_index{jj}),0,2)/sqrt(length(trial_index{jj})));
            Ulim(:,jj) = mnToPlot(:,jj)+(nanstd(elData(:,trial_index{jj}),0,2)/sqrt(length(trial_index{jj})));  
        end

        % Smooth the data?
        if specs.smoothLevel > 0
            temp = smooth(mnToPlot, specs.smoothLevel); mnToPlot = reshape(temp, size(mnToPlot));
            temp = smooth(Llim, specs.smoothLevel); Llim = reshape(temp, size(Llim));
            temp = smooth(Ulim, specs.smoothLevel); Ulim = reshape(temp, size(Ulim));
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
        
        % Plot max
        switch specs.plot.showMax
            case 'yes'
                timeInx = find(trials.time>0 & trials.time<1);
                [~,x] = max(abs(mnToPlot(timeInx,:))); 
                for jj = 1:length(trial_index)
                    p2 = plot(trials.time(timeInx(x(jj))),mnToPlot(timeInx(x(jj)),jj), 'Marker', '.', 'MarkerSize', 50, 'Color', colors(jj,:), 'LineStyle', 'none');
                    p2.Annotation.LegendInformation.IconDisplayStyle = 'off';        
                end
        end
         
%         if specs.plot.showMax > 0
%             timeInx = find(trials.time>0 & trials.time<1);
%             tmpMnToPlot = mnToPlot;
%             for mm = 1:specs.plot.showMax
%                 [~,x] = max(abs(tmpMnToPlot(timeInx,:))); 
%                 for jj = 1:length(trial_index)
%                     p2 = plot(trials.time(timeInx(x(jj))),mnToPlot(timeInx(x(jj)),jj), 'Marker', '.', 'MarkerSize', 50, 'Color', colors(jj,:), 'LineStyle', 'none');
%                     p2.Annotation.LegendInformation.IconDisplayStyle = 'off';        
%                 end
%                 tmpMnToPlot(timeInx(x),:) = nan;
%             end
%         end       
        
        % Check if these electrodes have matches with visual atlases, if so, add
        % that to the plot title
        isTableCol = @(t, thisCol) ismember(thisCol, t.Properties.VariableNames);
        if iscell(whichElectrodes)
            electrodeName = whichElectrodes{ii};
        else
            electrodeName = whichElectrodes;
        end
        viselec_name = [];
        if isfield(trials, 'viselec')
            for atlas = {'wang2015_atlas','benson14_varea', 'wang15_mplbl'}
                if ~isempty(trials.viselec)
                    if ~isempty(trials.viselec.(atlas{:}))
                        %viselec = contains(trials.viselec.(atlas{:}).elec_labels, electrodeName);
                        viselec = find(strcmp(electrodeName,trials.viselec.(atlas{:}).elec_labels));
                        if any(viselec)
                            %viselec_name = [viselec_name atlas{:}(1:8) ':' trials.viselec.(atlas{:}).area_labels{viselec} ' '];
                            atlasstr = strsplit(atlas{:},'_');
                            viselec_name = [viselec_name atlasstr{1}(1) ':' trials.viselec.(atlas{:}).area_labels{viselec} ' '];
                            switch atlas{:}
                                case 'benson14_varea'
                                    switch specs.plot.addEccToTitle
                                        case 'yes'
                                            viselec_name = [viselec_name sprintf('[ecc = %0.1f] ', trials.viselec.benson14_varea.node_eccen(viselec))];
                                             %viselec_name = [viselec_name ' ecc = ' num2str(trials.viselec.benson14_varea.node_eccen(viselec)) ' '];
                                     end
                            end                         
                        end
                    end
                end
            end
        else
            for atlas = {'wangarea','bensonarea'}
                if isTableCol(trials.channels,atlas{:})
                    viselec = find(strcmp(electrodeName,trials.channels.name));
                    visarea = trials.channels.(atlas{:}){viselec};
                    if ~isempty(visarea)&&~strcmp(visarea,'none')
                        atlasstr = strsplit(atlas{:},'_');
                        viselec_name = [viselec_name atlasstr{1}(1) ':' visarea ' '];
                        switch atlas{:}
                            case 'bensonarea'
                                switch specs.plot.addEccToTitle
                                    case 'yes'
                                        if isTableCol(trials.channels,'bensoneccen')
                                            viselec_name = [viselec_name sprintf('[ecc = %0.1f] ', trials.channels.bensoneccen{viselec})];
                                        end
                                end
                        end                         
                    end
                end
            end
        end
        
        % Set plot title
        plotTitle = [electrodeName ' ' viselec_name];
        out.titles{ii} = plotTitle;
        title(plotTitle);

        % Set y-axis limits
        if isempty(specs.plot.YLim)
            lim = [min(mnToPlot(:)) max(mnToPlot(:))];
            ylim = [lim(1)-(0.2*lim(1)*sign(lim(1))) lim(2)+(0.2*lim(2)*sign(lim(2)))];
        else
            ylim = specs.plot.YLim;
        end
        if any(isnan(ylim)), ylim = [-1 1]; end
        set(gca, 'YScale', specs.plot.YScale,'YLim', ylim);
        
        % Add legend
        if ii == subplotWithLegendIndex 
            switch specs.collapseTrialTypes
                    case 'no'
                        legend(trialType);
                    case 'yes'
                        legend([trialType{:}])
            end
        end
        
        % Add stim onset and zero lines
        l1 = line([0 0], ylim,'LineStyle', ':', 'Color', 'k');
        set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
        l2 = line([trials.time(1) trials.time(end)], [0 0],'LineStyle', ':', 'Color', 'k');
        set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
        
        % Save plotted lines in an output struct
        out.(thisDataType).(whichElectrodes{ii}).mn = double(mnToPlot');
        out.(thisDataType).(whichElectrodes{ii}).se = double((Ulim-mnToPlot)');
        set(gca, 'XScale', specs.plot.XScale, 'XLim', specs.plot.XLim);
        %set(gca, 'XLim', [-0.5 3]);
        
        % Set font size
        set(gca, 'FontSize', specs.plot.fontSize);
        
    end
    set(gcf, 'Position', [150 100 2000 1250]);
    %set(gcf, 'Position', [150 100 750 625]);
	out.time = trials.time;
    out.smooth = specs.smoothLevel;
    out.colors = colors;
    out.subplotdim = [nRow nCol];   
end

