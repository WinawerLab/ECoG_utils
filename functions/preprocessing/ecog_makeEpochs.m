function [epochs, t] = ecog_makeEpochs(raw_ts, onsets, epoch_time, fs)
% Slice time series matrix (channels x samples) into 3D epoched array
% (channels x samples x epochs) based on trigger times.
%
% [epochs, t] = ecog_makeEpochs(data, onsets, epoch_time, fs)
%
% Inputs
%   raw_ts:       time series matrix (channels x samples)
%   onsets:       a vector of stimulus onsets (either samples or seconds)
%   epoch_times:  a 2 vector of start and end time of the epochs 
%                   (in seconds) relative to the trial onset
%   fs:           sampling rate (Hz)
%
% Outputs
%   ts:           3D array containing epoched time series (samples by
%                   epoch x channel)

%% Parameters

if max(mod(onsets,1)) == 0
    % onsets are in samples
    onset_samples = onsets;
else
    % onsets are in seconds
    onset_samples = onsets*fs;
end
    
epoch_samples = round(epoch_time * fs); %epoch length in samples

inds = bsxfun(@plus,onset_samples,(epoch_samples(1):epoch_samples(2)-1));

ts       = raw_ts(:,inds); 
epochs   = reshape(ts, size(raw_ts,1),size(inds,1), size(inds,2));
epochs   = permute(epochs, [1 3 2]);

t = epoch_samples(1):epoch_samples(2)-1;
t = t/fs;

end
