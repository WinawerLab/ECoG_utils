function [ieeg_json] = bidsconvert_addChannelCounts(ieeg_json, channel_table);

% Update the ieeg_json file with the appropriate numbers from the
% channel_table

ieeg_json.ECOGChannelCount = length(find(contains(channel_table.type,'ecog')));
ieeg_json.SEEGChannelCount = length(find(contains(channel_table.type,'seeg')));
ieeg_json.EEGChannelCount = length(find(contains(channel_table.type,'eeg'))); 
ieeg_json.EOGChannelCount = length(find(contains(channel_table.type,'eog'))); 
ieeg_json.ECGChannelCount = length(find(contains(channel_table.type,'ecg'))); 
ieeg_json.EMGChannelCount = length(find(contains(channel_table.type,'emg'))); 
ieeg_json.MiscChannelCount = length([find(contains(channel_table.type,'misc')); find(contains(channel_table.type,'other')); find(contains(channel_table.type,'unknown'))]); 
ieeg_json.TriggerChannelCount = length(find(contains(channel_table.type,'trig'))); 