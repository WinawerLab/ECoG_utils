function ecog_plotSingleTimeCourse(t, timecourse, CI, lineSpec, plotTitle, yLabel, yLim)

% timecourse time x trials
% CI = timecourse * 2 [llim ulim]

if nargin < 7 || isempty(yLim)
    yLim = [];
end

if nargin < 6 || isempty(yLabel)
    yLabel= 'Broadband response';
end
    
if nargin < 5 || isempty(plotTitle)
    plotTitle= ' ';
end

if nargin < 4 || isempty(lineSpec)
    lineSpec = 'k';
end

if nargin < 3 || isempty(CI)
    CI = [];
end

if size(CI,2) ~= 2
    CI = CI';
end
 % Plot standard errors
if ~isempty(CI)
    h = ciplot(CI(:,1),CI(:,2),t,lineSpec, 0.25);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end

hold on;

% Plot means
if ischar(lineSpec)
    plot(t, timecourse, lineSpec, 'LineWidth',2);
else
    plot(t, timecourse, 'Color', lineSpec, 'LineWidth',2);
end

box off
axis tight

if isempty(yLim)
    yLim = get(gca, 'YLim');
else
    ylim(yLim);
end

% Add stim onset and zero lines
l1 = line([0 0], yLim,'LineStyle', ':', 'Color', 'k');
set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
l2 = line([t(1) t(end)], [0 0],'LineStyle', ':', 'Color', 'k');
set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 

% Add labels
xlabel('Time (s)'); 
ylabel(yLabel);
title(plotTitle);

set(gca, 'fontsize', 12);
 


end