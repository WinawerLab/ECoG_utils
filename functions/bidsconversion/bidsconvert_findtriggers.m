function [trigger_onsets,t] = bidsconvert_findtriggers(data, hdr, triggerChannel, makePlot, triggerIndex, peakOpts)

if nargin < 6
    peakOpts.minPeakHeight = 0.5;
    peakOpts.minPeakDistance = 0.5;
end

if nargin < 5 
    triggerIndex = [];
end

if nargin < 4 
    makePlot = 0;
end

% Get trigger time points from data file
triggers = data(triggerChannel,:);
triggers = triggers / max(triggers);
t = ((0:hdr.nSamples-1)/hdr.Fs); 

fprintf('[%s] Reading triggers from data \n', mfilename);
[~,trigger_onsets] = findpeaks(triggers, hdr.Fs, 'MinPeakHeight', peakOpts.minPeakHeight, 'MinPeakDistance', peakOpts.minPeakDistance);

if makePlot 
    figure('Name', 'Found triggers'); hold on;
    plot(t,triggers);
    if ~isempty(triggerIndex)
        plot(trigger_onsets(triggerIndex), ones(length(trigger_onsets(triggerIndex)),1),'r.','MarkerSize', 25, 'LineStyle','none');
    else
        plot(trigger_onsets, ones(length(trigger_onsets),1),'r.','MarkerSize', 25, 'LineStyle','none');
    end
    %stem(trigger_onsets,ones(length(trigger_onsets),1)*(max(data(triggerChannel,:))),'r');
    legend({'trigger data', 'trigger onsets'}); xlabel('time (s)'); ylabel('amplitude');
    title('found triggers');
    xlabel('Time (s)'); set(gca, 'FontSize', 12);
end
fprintf('[%s] Found %d triggers in total \n', mfilename, length(trigger_onsets));

end