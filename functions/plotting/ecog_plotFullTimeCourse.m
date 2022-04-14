function [out] = ecog_plotFullTimeCourse(data,dataType,chanstoplot,smoothLevel)

if nargin < 4 || isempty(smoothLevel)
    smoothLevel = 0;
end

switch data.channels.type{chanstoplot}
    case 'trig'
        data.(dataType)(chanstoplot(1),:) = data.(dataType)(chanstoplot(1),:)/max(data.(dataType)(chanstoplot(1),:));
        yval_onsets = max(data.(dataType)(chanstoplot(1),:));
    otherwise
        yval_onsets = mode(data.(dataType)(chanstoplot(:),:));
end

% Smooth & plot
figure; hold on
for ii = 1:length(chanstoplot)
    if smoothLevel > 0
        plot(data.time,smooth(data.(dataType)(chanstoplot(ii),:),smoothLevel), 'LineWidth', 2);
    else
        plot(data.time,data.(dataType)(chanstoplot(ii),:), 'LineWidth', 2);
    end
end

% Plot event onsets on top
plot(data.events.onset, ones(length(data.events.onset),1)*yval_onsets,'k.','MarkerSize', 25, 'LineStyle','none');
legend(data.hdr.label(chanstoplot));
xlabel('time (s)');
switch dataType
    case 'broadband'
        ylabel(['broadband power ( ' num2str(min(data.bb_bands(:))) '-' num2str(max(data.bb_bands(:))) ' Hz)']);
    case 'car_refer'
        ylabel('voltage (rereferenced to common average)');
    otherwise 
        ylabel('voltage');
end

l = line([data.time(1) max(data.time)], [0 0],'LineStyle', '-', 'Color', 'k');
set(gca, 'XLim', [data.time(1) data.time(end)]);
l.Annotation.LegendInformation.IconDisplayStyle = 'off';

end