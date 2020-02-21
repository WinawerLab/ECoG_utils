
function [out] = ecog_plotSpectra(spectra, whichElectrodes, trialType, taskType, specs)


if ~isfield(specs, 'plot') || isempty(specs.plot)
    specs.plot = [];
end

if ~isfield(specs.plot, 'colorMap') || isempty(specs.plot.colorMap), specs.plot.colorMap = 'copper'; end
if ~isfield(specs.plot, 'nSubPlots') || isempty(specs.plot.nSubPlots), specs.plot.nSubPlots = []; end
if ~isfield(specs.plot, 'addEccToTitle') || isempty(specs.plot.addEccToTitle), specs.plot.addEccToTitle = 'no'; end
if ~isfield(specs.plot, 'fontSize') || isempty(specs.plot.fontSize), specs.plot.fontSize = 12; end
if ~isfield(specs.plot, 'XLim') || isempty(specs.plot.XLim), specs.plot.XLim = [1 200];end
if ~isfield(specs.plot, 'YLim'), specs.plot.YLim = [];end

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
            trial_index{ii} = find(contains(spectra.events.trial_name, trialType{ii}) & contains(spectra.events.task_name, taskType));
        end
    else
        for ii = 1:length(trialType)
            trial_index{ii} = find(contains(spectra.events.trial_name, trialType{ii}));
        end
    end
end

% Plot settings
cmap = eval(specs.plot.colorMap);
nCond = length(trial_index);
colors = cmap(1:round((length(cmap)/nCond)):length(cmap),:);
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

figure('Name', [trialType{~contains(trialType, 'BLANK')}]); 

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
            mnToPlot(:,jj) = squeeze(mean(elData(:,trial_index{jj}),2));
            Llim(:,jj) = mnToPlot(:,jj)-(std(elData(:,trial_index{jj}),0,2)/sqrt(length(trial_index{jj})));
            Ulim(:,jj) = mnToPlot(:,jj)+(std(elData(:,trial_index{jj}),0,2)/sqrt(length(trial_index{jj})));        
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
        set(gca, 'Yscale', 'log','YLim', ylim);
        %set(gca,'YLim', ylim);
        % Set x-axis limits
        set(gca, 'XLim', specs.plot.XLim);
        
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
set(gcf, 'Position', [150 100 2000 1250]);

out.f           = spectra.f;
out.colors      = colors;
out.subplotdim  = [nRow nCol];   

end