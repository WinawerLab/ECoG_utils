function [epochs, t] = ecog_makeEpochs(raw_ts, onsets, epoch_time, fs)
% Slice time series matrix (channels x samples) into 3D epoched array
% (samples x epochs x channels) based on trigger times.
%
% [epochs, t] = ecog_makeEpochs(data, onsets, epoch_time, fs)
%
% Inputs:
%   raw_ts:       time series matrix (channels x samples)
%   onsets:       a vector of stimulus onsets (in samples or seconds)
%   epoch_times:  a 2 vector of start and end time of the epochs 
%                   (in seconds) relative to the trial onset
%   fs:           sampling rate (Hz)
%
% Outputs:
%   epochs:       3D array containing epoched time series (samples x
%                   events x channels)                   
%   t:            1D array of epoch time points relative to stimulus onset 
%                   e.g. -0.5 to 1 (in seconds)


%% Parameters

% Determine if the onsets are provided in samples or seconds
if max(mod(onsets,1)) == 0
    % onsets are in samples
    onset_samples = onsets;
else
    % onsets are in seconds
    onset_samples = round(onsets*fs);
end
    
epoch_samples = round(epoch_time * fs); %epoch length in samples

%% Create epochs

inds = bsxfun(@plus,onset_samples,(epoch_samples(1):epoch_samples(2)-1));

ts       = raw_ts(:,inds); 
epochs   = reshape(ts, size(raw_ts,1),size(inds,1), size(inds,2));
epochs   = permute(epochs, [3 2 1]);

t = epoch_samples(1):epoch_samples(2)-1;
t = t/fs;
t = t';

end
