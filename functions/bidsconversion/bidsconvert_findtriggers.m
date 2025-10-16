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

%% trigger selection for tactile task

[~,trigger_onsets, widths] = findpeaks(triggers, hdr.Fs,...
    'MinPeakHeight',peakOpts.minPeakHeight,...
    'MinPeakProminence',peakOpts.minPeakHeight,...
    'MinPeakDistance', peakOpts.minPeakDistance);
valid_index = widths < 0.1;
trigger_onsets = trigger_onsets(valid_index);

% Segment triggers into blocks using the gap threshold (e.g., 5 seconds)
dt = diff(trigger_onsets);                   % Inter-trigger intervals
block_break_indices = find(dt > 5);            % Find large gaps indicating block boundaries
numBlocks = length(block_break_indices) + 1;   % Number of blocks is one more than the number of breaks

% Split trigger_onsets into blocks
blocks = cell(numBlocks, 1);
startIdx = 1;
for i = 1:length(block_break_indices)
    endIdx = block_break_indices(i);
    blocks{i} = trigger_onsets(startIdx:endIdx);
    startIdx = endIdx + 1;
end
blocks{numBlocks} = trigger_onsets(startIdx:end);

% Display block counts before correction
fprintf('Trigger counts per block BEFORE correction:\n');
for i = 1:numBlocks
    fprintf('Block %d: %d triggers\n', i, numel(blocks{i}));
end

% % For blocks with 145 triggers, remove the first trigger
% for i = 1:numBlocks
%     if numel(blocks{i}) == 145
%         blocks{i}(1) = [];  % Remove the first trigger
%         fprintf('Block %d had 145 triggers: first trigger removed.\n', i);
%     end
% end

% Display block counts after correction
fprintf('\nTrigger counts per block AFTER correction:\n');
for i = 1:numBlocks
    fprintf('Block %d: %d triggers\n', i, numel(blocks{i}));
end

% Combine the corrected blocks back into a single ordered vector
corrected_trigger_onsets = horzcat(blocks{:});
% Final total trigger count (expecting 144 * number of blocks)
fprintf('\nTotal triggers after correction: %d\n', numel(corrected_trigger_onsets));
trigger_onsets = corrected_trigger_onsets;

%% trigger selection for other tasks
% [~,trigger_onsets] = findpeaks(triggers, hdr.Fs,'MinPeakHeight',peakOpts.minPeakHeight,'MinPeakDistance', peakOpts.minPeakDistance);

%% check trigger
if makePlot
    figure('Name', 'Found triggers'); hold on;
    plot(t,triggers);
    if ~isempty(triggerIndex)
        plot(trigger_onsets(triggerIndex), ones(length(trigger_onsets(triggerIndex)),1),'r.','MarkerSize', 25, 'LineStyle','none');
    else
        plot(trigger_onsets, ones(length(trigger_onsets),1),'r.','MarkerSize', 25, 'LineStyle','none');
    end
    % stem(trigger_onsets,ones(length(trigger_onsets),1)*(max(data(triggerChannel,:))),'r');
    legend({'trigger data', 'trigger onsets'}); xlabel('time (s)'); ylabel('amplitude');
    title('found triggers');
    xlabel('Time (s)'); set(gca, 'FontSize', 12);
end
fprintf('[%s] Found %d triggers in total \n', mfilename, length(trigger_onsets));

end