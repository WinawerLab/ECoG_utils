function [trials_out,figlist] = ecog_plotGridTimecourses(trials, whichGrid, whichTrials, specs)

% [trials_out, p] = ecog_plotGridTimecourses(trials, whichGrid, whichTrials, specs)
% Plot in the entire grid.
%   trials.evoked / trials.broadband = chs x times x trials
% 
% See also ecog_plotTimecourses

% Dependency: plotGridSimpleCommon, SetDefault

% 20190725 Yuasa
% 20200226 Yuasa: output figure handles
% 20220222 Yuasa: use plotGridSimpleCommon

%% Set options
narginchk(3,inf);
SetDefault('specs.dataTypes',{'broadband', 'evoked'}, 0);
SetDefault('specs.plot.XScale','linear', 0);
SetDefault('specs.plot.YScale','linear', 0);
SetDefault('specs.plot.nSubPlots', [], 1);
SetDefault('specs.plot.RotGrid', false, 0);
SetDefault('specs.plot.showlegend','Outside', 0);    % Outside, Inside, Last, No

%%
%%% run plotGridSimpleCommon
channels = trials.channels;
plotGridSimpleCommon;
trials.channels = channels;

%-- Plot figures
for ee = 1:length(inx)
    
    [trials_out{ee}] = ecog_plotTimecourses(trials, gridList(inx{ee}), whichTrials, specs);
    
end

%-- Set axis scale
figlist = get(groot,'Children');
figlist = flipud(figlist(1:length(inx)*length(specs.dataTypes)))';
for ifig=figlist
    axlist = findobj(ifig,'Type','Axes')';
    for iax = axlist
        if ~isempty(iax.Title.String)
            iax.XScale = specs.plot.XScale;
            iax.YScale = specs.plot.YScale;
            if strcmp(iax.XScale,'log')
                if numel(iax.XTick) < 2
                    curtick = [log10(iax.XTick) floor(log10(max(trials.time)))];
                    if diff(curtick) == 0
                        iax.XTick = 10.^([-2:0] + curtick(1));
                    else
                        iax.XTick = 10.^([curtick(1)-diff(curtick) curtick]);
                    end
                elseif numel(iax.XTick) < 3
                    curtick = log10(iax.XTick);
                    iax.XTick = 10.^([curtick(1), mean(curtick), curtick(end)]);
                end
                iax.XTickLabel = cellstr(num2str(iax.XTick'));
            end
        end
    end  
end

%-- Hide or reorder legend
shiftlefend(figlist,specs.plot.showlegend);

end

%% Sub function
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
if ~isfield(specs.plot, 'XLim') || isempty(specs.plot.XLim), specs.plot.XLim = [];end
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
            elData = trials.(thisDataType)(el_index(ii),:,:);
            elData = permute(elData,[2:ndims(elData),1]);   % squeeze works wrong when data was averaged
        else
            subplotWithLegendIndex = subplotWithLegendIndex+1;
            set(gca, 'FontSize', specs.plot.fontSize);
            continue
        end
        
        % Baseline correction: 
        switch specs.baselineType
            case 'all'
                baseline_index = find(~contains(trials.events.trial_name, {'PRF', 'BLANK'}));
                if isempty(baseline_index), baseline_index = true(height(trials.events),1); end
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
            Llim(:,jj) = mnToPlot(:,jj)-(nanstd(elData(:,trial_index{jj}),0,2)./sqrt(sum(~isnan(elData(:,trial_index{jj})),2)));
            Ulim(:,jj) = mnToPlot(:,jj)+(nanstd(elData(:,trial_index{jj}),0,2)./sqrt(sum(~isnan(elData(:,trial_index{jj})),2)));
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
        if isempty(specs.plot.XLim)
            lim = [min(mnToPlot(:)) max(mnToPlot(:))];
            ylim = [lim(1)-(0.2*lim(1)*sign(lim(1))) lim(2)+(0.2*lim(2)*sign(lim(2)))];
        else
            xlim = specs.plot.XLim;
        end
        out.(thisDataType).(whichElectrodes{ii}).mn = double(mnToPlot');
        out.(thisDataType).(whichElectrodes{ii}).se = double((Ulim-mnToPlot)');
        if isempty(specs.plot.XLim)
            set(gca, 'XScale', specs.plot.XScale);
        else
            set(gca, 'XScale', specs.plot.XScale, 'XLim', specs.plot.XLim);
        end
        
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
end

function shiftlefend(figlist,showlegend)

switch lower(showlegend)
    case 'no'
        set(findobj(figlist,'Type','Legend'),'Visible','off');  % Hide
    case {'outside','last'}
        for ifig=figlist
            %%% get axis size
            hax   = findobj(ifig,'Type','Axes'); 
            frstaxopos = hax(end).OuterPosition;
            frstaxpos  = hax(end).Position;
            lastaxopos = hax(1).OuterPosition;
            lastaxpos  = hax(1).Position;
            
            %%% get index of legend
            hleg  = findobj(ifig,'Type','legend');
            hleg  = hleg(1);
            
            %%% move legend at outer position
            legpos = hleg.Position;
            if strcmpi(showlegend,'last')
                hleg.Position = ...
                    [sum(lastaxopos([1,3]))+lastaxopos(1)./10,...
                     lastaxpos(2) + lastaxpos(4)./2 - legpos(4)./2,...
                     legpos(3:4)];
            else
                hleg.Position = ...
                    [frstaxpos(1)+(frstaxpos(3)-legpos(3))/2, frstaxopos(2)+frstaxopos(4), legpos(3:4)];
            end
            
            %%% move legend and axes next to legend at the top of axes list
            legidx = find(ifig.Children==hleg);
            uistack([hleg,ifig.Children(legidx+1)],'top');
        end
end
end
