function ecog_plotMultipleTimeCourses(t, timecourses, CI, colors, plotTitle, yLabel, yLim, xLabel, xLim)

% Plots a response timecourse 
%
% ecog_plotSingleTimeCourse(t, timecourse, [CI], [lineSpec], [plotTitle], [yLabel], [yLim], [xLabel], [xLim])
%
% INPUTS
%   t:              xvalues (time)
%   timecourses:    yvalues (e.g. time-varying broadband or ERP) time x
%                   timecourses
%   CI:             (optional) timecourses reflecting confidence intervals 
%                   should be size of timecourses * 2 [lowerlim upperlim]
%                       default: no CI
%   colors:         (optional) color label (e.g., 'r') or RGB triplet
%                       default: 'k'
%   plotTitle:      (optional) title of plot
%   yLabel:         (optional) y-axis label
%   yLim:           (optional) min and max of y-axis
%   xLabel:         (optional) x-axis label
%   xLim:           (optional) min and max of x-axis

[~, nTs] = size(timecourses);

if ~exist('CI', 'var'), CI = []; end
if ~exist('colors', 'var') || isempty(colors), colors = parula(nTs); end
if ~exist('plotTitle', 'var'), plotTitle= []; end
if ~exist('yLabel', 'var'), yLabel = []; end
if ~exist('yLim', 'var'), yLim = []; end
if ~exist('xLabel', 'var'), xLabel = []; end
if ~exist('xLim', 'var'), xLim = []; end

% Check if CI is formatted correctly
if ~isempty(CI) && size(CI,3) ~= 2
    fprintf('[%s] CI not formatted correctly for plotting\n', mfilename);
    return
end


% Plot CI first

for ii = 1:nTs
    
    if ~isempty(CI)
        h = ciplot(CI(:,ii,1),CI(:,ii,2),t, colors, 0.25);
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end

hold on;

for ii = 1:nTs

    % Plot time course
    if iscell(colors)
        plot(t, timecourses(:,ii), colors{ii}, 'LineWidth',2);
    else
        plot(t, timecourses(:,ii), 'Color', colors(ii,:), 'LineWidth',2);
    end
end

box off
axis tight

% Set axes
if ~isempty(xLim), xlim(xLim); end

if isempty(yLim)
    yLim = get(gca, 'YLim');
else
    ylim(yLim);
end

% add stim onset and zero lines
%l1 = line([0 0], yLim,'LineStyle', ':', 'Color', 'k');
l1 = line([0 0], [-10^3 10^3],'LineStyle', ':', 'Color', 'k');
set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
l2 = line([t(1) t(end)], [0 0],'LineStyle', ':', 'Color', 'k');
set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 

% Add labels
if ~isempty(xLabel), xlabel(xLabel); end
if ~isempty(yLabel), ylabel(yLabel); end
if ~isempty(plotTitle), title(plotTitle); end
set(gca, 'yLim',yLim);
set(gca, 'fontsize',14);
 
end