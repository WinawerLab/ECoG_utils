function ecog_plotSingleTimeCourse(t, timecourse, CI, lineSpec, plotTitle, yLabel, yLim, xLabel, xLim)

% Plots a response timecourse 
%
% ecog_plotSingleTimeCourse(t, timecourse, [CI], [lineSpec], [plotTitle], [yLabel], [yLim], [xLabel], [xLim])
%
% INPUTS
%   t:              xvalues (time)
%   timecourse:     yvalues (e.g. time-varying broadband or ERP)
%   CI:             (optional) timecourses reflecting confidence intervals 
%                   should be size of timecourse * 2 [lowerlim upperlim]
%                       default: no CI
%   linespec:       (optional) color label (e.g., 'r') or RGB triplet
%                       default: 'k'
%   plotTitle:      (optional) title of plot
%   yLabel:         (optional) y-axis label
%   yLim:           (optional) min and max of y-axis
%   xLabel:         (optional) x-axis label
%   xLim:           (optional) min and max of x-axis


if ~exist('CI', 'var'), CI = []; end
if ~exist('lineSpec', 'var') || isempty(lineSpec), lineSpec = 'k'; end
if ~exist('plotTitle', 'var'), plotTitle= []; end
if ~exist('yLabel', 'var'), yLabel = []; end
if ~exist('yLim', 'var'), yLim = []; end
if ~exist('xLabel', 'var'), xLabel = []; end
if ~exist('xLim', 'var'), xLim = []; end

% Check if CI is formatted correctly
if size(CI,2) ~= 2
    CI = CI';
end

% Plot CI first
if ~isempty(CI)
    h = ciplot(CI(:,1),CI(:,2),t,lineSpec, 0.25);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end

hold on;

% Plot time course
if ischar(lineSpec)
    plot(t, timecourse, lineSpec, 'LineWidth',2);
else
    plot(t, timecourse, 'Color', lineSpec, 'LineWidth',2);
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

% Plot stim onset and zero lines
h = get(gcf);
if isempty(h.CurrentObject)
    l1 = line([0 0], yLim,'LineStyle', ':', 'Color', 'k');
    %l1 = line([0 0], [-10^4 10^4],'LineStyle', ':', 'Color', 'k');
    set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
    l2 = line([t(1) t(end)], [0 0],'LineStyle', ':', 'Color', 'k');
    set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
end

% Add labels
if ~isempty(xLabel), xlabel(xLabel); end
if ~isempty(yLabel), ylabel(yLabel); end
if ~isempty(plotTitle), title(plotTitle); end
set(gca, 'fontsize',14);
 

end