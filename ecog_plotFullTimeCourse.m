function [out] = ecog_plotFullTimeCourse(timecourse,chanstoplot,events)

if nargin < 3 || isempty(events)
    events.onset = [];
end
    
% Smooth & plot
figure; hold on
for ii = 1:length(chanstoplot)
    plot(data.time{1},smooth(timecourse(chanstoplot(ii),:),128), 'LineWidth', 2);
end

% and plot event onsets on top
plot(events.onset, zeros(length(events.onset),1),'k.','MarkerSize', 25, 'LineStyle','none');
legend(data.label(chanstoplot));
xlabel('time (s)');
ylabel(['broadband power ( ' num2str(min(min(bands))) '-' numstr(max(max(bands))) ' Hz)']);
%ylabel('evoked');
l = line([data.time{1}(1) max(data.time{1})], [0 0],'LineStyle', '-', 'Color', 'k');
l.Annotation.LegendInformation.IconDisplayStyle = 'off';

end