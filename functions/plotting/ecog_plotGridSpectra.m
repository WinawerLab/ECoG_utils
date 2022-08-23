function [trials_out,figlist] = ecog_plotGridSpectra(spectra, whichGrid, whichTrials, whichTasks, specs)

% [trials_out, p] = ecog_plotGridSpectra(trials, whichGrid, whichTrials, whichTasks, specs)
% Plot in the entire grid.
% 
% See also ecog_plotSpectra

% Dependency: plotGridSimpleCommon, SetDefault

% 20200306 Yuasa
% 20220222 Yuasa: use plotGridSimpleCommon

%% Set options
narginchk(3,inf);
SetDefault('specs.plot.XScale','linear', 0);
SetDefault('specs.plot.YScale','log', 0);
SetDefault('specs.plot.nSubPlots', [], 1);
SetDefault('specs.plot.RotGrid', false, 0);
SetDefault('specs.plot.showlegend','Outside', 0);    % Outside, Inside, Last, No

%%
%%% run plotGridSimpleCommon
channels = spectra.channels;
plotGridSimpleCommon;
spectra.channels = channels;

%-- Plot figures
for ee = 1:length(inx)
    
    [trials_out{ee}] = ecog_plotSpectra(spectra, gridList(inx{ee}), whichTrials, whichTasks, specs);
    
end

%-- Set axis scale
figlist = get(groot,'Children');
figlist = flipud(figlist(1:length(inx)))';
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
                iax.XTickLabel = cellstr(num2str(iax.XTick',['%.' int2str(log10(diff(iax.XLim))-2) 'f']));
            end
        end
    end  
end

%-- Hide or reorder legend
shiftlefend(figlist,specs.plot.showlegend);

end

%% Sub function
function [out] = ecog_plotSpectra(spectra, whichElectrodes, trialType, taskType, specs)


if ~isfield(specs, 'plot') || isempty(specs.plot)
    specs.plot = [];
end

if ~isfield(specs, 'trialIntpr') || isempty(specs.trialIntpr)
    specs.trialIntpr = @contains;
end
if ~isfield(specs, 'taksIntpr') || isempty(specs.taksIntpr)
    specs.taksIntpr = @contains;
end

if ~isfield(specs.plot, 'colorMap') || isempty(specs.plot.colorMap), specs.plot.colorMap = 'copper'; end
if ~isfield(specs.plot, 'nSubPlots') || isempty(specs.plot.nSubPlots), specs.plot.nSubPlots = []; end
if ~isfield(specs.plot, 'addEccToTitle') || isempty(specs.plot.addEccToTitle), specs.plot.addEccToTitle = 'no'; end
if ~isfield(specs.plot, 'fontSize') || isempty(specs.plot.fontSize), specs.plot.fontSize = 12; end
if ~isfield(specs.plot, 'XLim') || isempty(specs.plot.XLim), specs.plot.XLim = [1 200];end
if ~isfield(specs.plot, 'YLim'), specs.plot.YLim = [];end
if ~isfield(specs.plot, 'XScale') || isempty(specs.plot.XScale), specs.plot.XScale = 'linear';end
if ~isfield(specs.plot, 'YScale') || isempty(specs.plot.YScale), specs.plot.YScale = 'log';end

out = struct;

% Find electrodes in data
elInx = ecog_matchChannels(whichElectrodes, spectra);
if isempty(elInx)
    fprintf('[%s] Did not find any matching channels, not plotting \n', mfilename);  
    return
end

% Find trials in data
trial_index = [];
if isempty(trialType)
    disp('Please specify a trial type')
    return
else
    if ~isempty(taskType)
        for ii = 1:length(trialType)
            trial_index{ii} = find(specs.trialIntpr(spectra.events.trial_name, trialType{ii}) & specs.taksIntpr(spectra.events.task_name, taskType));
        end
    else
        for ii = 1:length(trialType)
            trial_index{ii} = find(specs.trialIntpr(spectra.events.trial_name, trialType{ii}));
        end
    end
end

% Plot settings
cmap = eval(specs.plot.colorMap);
nCond = length(trial_index);
colors = cmap(round(1:(length(cmap)/nCond):length(cmap)),:);
%colors = sortrows(colors,'descend');
if max(contains(trialType, 'BLANK') >0)
    nCond = nCond-1; 
    colors(1,:) = [];
end

out.elecs = whichElectrodes;
out.titles = cell(size(out.elecs));

% Decide how many subplots are needed
if ~isempty(specs.plot.nSubPlots)
    nRow = specs.plot.nSubPlots(1);
    nCol = specs.plot.nSubPlots(2);
else
    nPlot = length(elInx);
    nRow = ceil(sqrt(nPlot));
    nCol = ceil(sqrt(nPlot));
    if nPlot <= (nRow*nCol)-nCol
        nRow = nRow-1;
    end
end

%figureName = strsplit(trialType{1},'-');
%figure('Name', [figureName{1} ' ' thisDataType]); 

figure('Name', [trialType{~contains(trialType, 'BLANK')}],'Position', [150 100 2000 1250]); 

hasLegend = 0;
for ee = 1:length(elInx)
    
    subplot(nRow,nCol,ee); hold on;
    
    % Set font size
    set(gca, 'FontSize', specs.plot.fontSize);
    
    if elInx(ee) > 0
        
        elData = squeeze(spectra.pwrspctrm(elInx(ee),:,:))';

        % Select subset of trials to plot
        clear mnToPlot Llim Ulim
        for jj = 1:length(trial_index)
            % Compute mean and standard error of the mean
            mnToPlot(:,jj) = squeeze(nanmean(elData(:,trial_index{jj}),2));
            Llim(:,jj) = mnToPlot(:,jj)-(nanstd(elData(:,trial_index{jj}),0,2)/sqrt(length(trial_index{jj})));
            Ulim(:,jj) = mnToPlot(:,jj)+(nanstd(elData(:,trial_index{jj}),0,2)/sqrt(length(trial_index{jj})));        
        end

%         % Smooth the data?
%         if specs.smoothLevel > 0
%             temp = smooth(mnToPlot, specs.smoothLevel); mnToPlot = reshape(temp, size(mnToPlot));
%             temp = smooth(Llim, specs.smoothLevel); Llim = reshape(temp, size(Llim));
%             temp = smooth(Ulim, specs.smoothLevel); Ulim = reshape(temp, size(Ulim));
%         end
        
        % Plot standard errors
        colorInx = 1;
        for jj = 1:length(trial_index)
            if strcmp(trialType{jj}, 'BLANK')
                colortoplot = [0 0 0];
            else
                colortoplot = colors(colorInx,:); colorInx = colorInx + 1;
            end
            h = ciplot(Llim(:,jj),Ulim(:,jj),spectra.f,colortoplot, 0.25);
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        end
        
        % Plot means
        colorInx = 1;
        for jj = 1:length(trial_index)
            if strcmp(trialType{jj}, 'BLANK')
                colortoplot = [0 0 0];
            else
                colortoplot = colors(colorInx,:); colorInx = colorInx + 1;
            end
            plot(spectra.f, mnToPlot(:,jj),'Color', colortoplot, 'LineWidth',2);
        end
              
        % Check if these electrodes have matches with visual atlases, if so, add
        % that to the plot title
        isTableCol = @(t, thisCol) ismember(thisCol, t.Properties.VariableNames);
        if iscell(whichElectrodes)
            electrodeName = whichElectrodes{ee};
        else
            electrodeName = whichElectrodes;
        end
        viselec_name = [];
        if isfield(spectra, 'viselec')
            for atlas = {'wang2015_atlas','benson14_varea', 'wang15_mplbl'}
                if ~isempty(spectra.viselec)
                    if ~isempty(spectra.viselec.(atlas{:}))
                        viselec = find(strcmp(electrodeName,spectra.viselec.(atlas{:}).elec_labels));
                        if any(viselec)
                            atlasstr = strsplit(atlas{:},'_');
                            viselec_name = [viselec_name atlasstr{1}(1) ':' spectra.viselec.(atlas{:}).area_labels{viselec} ' '];
                            switch atlas{:}
                                case 'benson14_varea'
                                    switch specs.plot.addEccToTitle
                                        case 'yes'
                                            viselec_name = [viselec_name sprintf('[ecc = %0.1f] ', spectra.viselec.benson14_varea.node_eccen(viselec))];
                                     end
                            end                         
                        end
                    end
                end
            end
        else
            for atlas = {'wangarea','bensonarea'}
                if isTableCol(spectra.channels,atlas{:})
                    viselec = find(strcmp(electrodeName,spectra.channels.name));
                    visarea = spectra.channels.(atlas{:}){viselec};
                    if ~isempty(visarea)&&~strcmp(visarea,'none')
                        atlasstr = strsplit(atlas{:},'_');
                        viselec_name = [viselec_name atlasstr{1}(1) ':' visarea ' '];
                        switch atlas{:}
                            case 'bensonarea'
                                switch specs.plot.addEccToTitle
                                    case 'yes'
                                        if isTableCol(spectra.channels,'bensoneccen')
                                            viselec_name = [viselec_name sprintf('[ecc = %0.1f] ', spectra.channels.bensoneccen{viselec})];
                                        end
                                end
                        end                         
                    end
                end
            end
        end
        
        % Set plot title
        plotTitle = [electrodeName ' ' viselec_name];
        out.titles{ee} = plotTitle;
        title(plotTitle);
        
        % Set y-axis limits
        if isempty(specs.plot.YLim)
            lim = [min(mnToPlot(:)) max(mnToPlot(:))];
            ylim = [lim(1)-(0.2*lim(1)*sign(lim(1))) lim(2)+(0.2*lim(2)*sign(lim(2)))];
        else
            ylim = specs.plot.YLim;
        end
        if any(isnan(ylim)), ylim = [0 1]; end
        set(gca, 'YScale', specs.plot.YScale,'YLim', ylim);
        % Set x-axis limits
        set(gca, 'XScale', specs.plot.XScale, 'XLim', specs.plot.XLim);
        if all(xticks>(specs.plot.XLim(1)+diff(specs.plot.XLim)/10))
            xticks([specs.plot.XLim(1), xticks]);
        end
        % Avoid exponential notation
        set(gca, 'XTickLabel', cellstr(num2str(get(gca,'XTick')',['%.' int2str(max(0,1-log10(diff(specs.plot.XLim)))) 'f'])));
        
        % Add legend
        if hasLegend == 0
            legend(trialType);
            xlabel('frequency');
            ylabel('power');
            hasLegend = 1;
        end  
        
        out.pwrspctrm.(whichElectrodes{ee}).mn = double(mnToPlot');
        out.pwrspctrm.(whichElectrodes{ee}).se = double((Ulim-mnToPlot)');
    end
end

out.f           = spectra.f;
out.colors      = colors;
out.subplotdim  = [nRow nCol];   

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
