function ecog_plotSingleSpectrum(f, spectrum, CI, lineSpec, yLabel, yLim, yScale, xLabel, xLim, xScale, plotTitle)

% Plots a response spectrum
%
% ecog_plotSingleSpectrum(f, spectrum, [CI], [lineSpec], [yLabel], [yLim], [yScale], [xLabel], [xLim], [xScale], [plotTitle])
%
% INPUTS
%   t:              xvalues (time)
%   timecourse:     yvalues (e.g. time-varying broadband or ERP)
%   CI:             (optional) timecourses reflecting confidence intervals 
%                   should be size of timecourse * 2 [lowerlim upperlim]
%                       default: no CI
%   linespec:       (optional) color label (e.g., 'r') or RGB triplet
%                       default: 'k'
%   yLabel:         (optional) y-axis label
%   yLim:           (optional) min and max of y-axis
%   yScale:         (optional) log or linear (default log)
%   xLabel:         (optional) x-axis label
%   xLim:           (optional) min and max of x-axis
%   xScale:         (optional) log or linear (default linear)
%   plotTitle:      (optional) title of plot

if ~exist('CI', 'var'), CI = []; end
if ~exist('lineSpec', 'var') || isempty(lineSpec), lineSpec = 'k'; end
if ~exist('yLabel', 'var'), yLabel = []; end
if ~exist('yLim', 'var'), yLim = []; end
if ~exist('yScale', 'var'), yScale = 'log'; end
if ~exist('xLabel', 'var'), xLabel = []; end
if ~exist('xLim', 'var'), xLim = []; end
if ~exist('xScale', 'var'), xScale = 'linear'; end
if ~exist('plotTitle', 'var'), plotTitle= []; end

% Check if CI is formatted correctly
if size(CI,2) ~= 2
    CI = CI';
end

% Plot standard errors first
if ~isempty(CI)
    h = ciplot(CI(:,1),CI(:,2),f,lineSpec, 0.25);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end

hold on;

% Plot time course
if ischar(lineSpec)
    plot(f, spectrum, lineSpec, 'LineWidth',2);
else
    plot(f, spectrum, 'Color', lineSpec, 'LineWidth',2);
end

box off
axis tight

% Set x-axis 
set(gca, 'XScale', xScale);
if ~isempty(xLim), xlim(xLim); end

%if all(xticks>(specs.plot.XLim(1)+diff(specs.plot.XLim)/10))
%    xticks([specs.plot.XLim(1), xticks]);
%end

% Set y-axis limits
set(gca, 'YScale', yScale);
if isempty(yLim)
    lim = get(gca, 'YLim');
    yLim = [lim(1)-(0.2*lim(1)*sign(lim(1))) lim(2)+(0.2*lim(2)*sign(lim(2)))];
end
ylim(yLim);
     

% Add labels
if ~isempty(xLabel), xlabel(xLabel); end
if ~isempty(yLabel), ylabel(yLabel); end
if ~isempty(plotTitle), title(plotTitle); end

set(gca, 'fontsize',10);
 


end