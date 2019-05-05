function [trigger_onsets] = bidsconvert_findtriggers(data, hdr, triggerChannel, makePlot)

if nargin < 4 
    makePlot = 0;
end

% Get trigger time points from data file
triggers = data(triggerChannel,:);
triggers = triggers / max(triggers);
t = ((0:hdr.nSamples-1)/hdr.Fs); 

fprintf('[%s] Reading triggers from data \n', mfilename);
[~,trigger_onsets] = findpeaks(triggers, hdr.Fs, 'MinPeakHeight', 0.8, 'MinPeakDistance', 0.5);

if makePlot 
    figure('Name', 'Found triggers'); hold on;
    plot(t,triggers);
    plot(trigger_onsets, ones(length(trigger_onsets),1),'r.','MarkerSize', 25, 'LineStyle','none');
    %stem(trigger_onsets,ones(length(trigger_onsets),1)*(max(data(triggerChannel,:))),'r');
    legend({'trigger data', 'trigger onsets'}); xlabel('time (s)'); ylabel('amplitude');
    title('found triggers');
    xlabel('Time (s)'); set(gca, 'FontSize', 12);
end
fprintf('[%s] Found %d triggers in total \n', mfilename, length(trigger_onsets));

end