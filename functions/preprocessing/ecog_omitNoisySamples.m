function [data,omits] = ecog_omitNoisySamples(data,srate,noise_freq,nharmonic,alpha,minsamples,gapsamples,isfill)

% data = ecog_omitNoisySamples(data, srate, noise_freq, [nharmonic, alpha, minsamples, gapsamples])
%   fills NaN values in time points where artifact noise is significantly large.
%
% [data, omits] = ecog_omitNoisySamples(data, srate, noise_freq, ...)
%   returns 'omits', a logical array indicating omitted samples.
%
% data = ecog_omitNoisySamples(..., 'fill')
%   fills missing values after omitting noisy samples.
%
% INPUTS:
%   data:        (time X channels)
%   srate:       Sampling rate (Hz)
%   noise_freq:  Artifact noise frequency (Hz) [default = 60]
%   nharmonic:   Nth harmonics to be notched [default = 3]
%   alpha:       Exclusion threshold [default = 0.003 (0.3%)]
%   minsamples:  Omit samples if noise appears for at least N samples [default = 3]
%   gapsamples:  Omit N samples that get caught in noise samples [default = 2]

% 20230816 Yuasa

%-- Set parameters
narginchk(2,inf);
if ~exist('noise_freq', 'var') || isempty(noise_freq), noise_freq = 60; end  % noise frequency
if ~exist('nharmonic', 'var') || isempty(nharmonic), nharmonic = 3; end      % Nth-harmonic
if ~exist('alpha', 'var') || isempty(alpha), alpha = 0.003; end              % [0,1]
if ~exist('minsamples', 'var') || isempty(minsamples), minsamples = 3; end   % Integar
if ~exist('gapsamples', 'var') || isempty(gapsamples), gapsamples = 2; end   % Integar
isfill = exist('isfill', 'var') && ischar(isfill) && strcmpi(isfill,'fill');  % 'fill' or empty

nharmonic  = round(nharmonic(1));
gapsamples = round(gapsamples(1));

%-- Expand notch frequencies to include all harmonics
noise_freq = reshape(noise_freq(:).*(1:nharmonic),[],1) + [-5 5];

%-- Compute power timecourse at noise frequnecy
fprintf('[%s] Omitting noisy samples ...\n', mfilename);
[noise_signal] = ecog_extractBroadband(data, srate, [], noise_freq);

%-- Fill NaN for each channel
omits   = false(size(data));
for ich = 1:size(noise_signal,2)
    inoise = noise_signal(:,ich);
    estdist = fitdist(inoise,'Loglogistic');
    noise_samples = inoise > estdist.icdf(1-alpha);
    %--- Ignore small gap 
    noise_clstr   = cell2mat(arrayfun(@(t)t.*ones(t,1),diff(find([1;diff(noise_samples);1])),'UniformOutput',false));
    noise_samples(noise_clstr.*~noise_samples <= gapsamples) = true;
    %--- Ignore short noise 
    noise_clstr   = cell2mat(arrayfun(@(t)t.*ones(t,1),diff(find([1;diff(noise_samples);1])),'UniformOutput',false));
    noise_samples(noise_clstr.*noise_samples < minsamples) = false;
    %-- Fill NaN
    data(noise_samples,ich) = nan;
    omits(:,ich)            = noise_samples;
end
omitrat = sum(omits,1)./size(omits,1);

%-- Filling missing
if isfill
    datsiz = size(data);
    padding = zeros([2,size(data,2:ndims(data))]);
    data = cat(1,padding,data,padding);
    data = fillmissing(data,'pchip',1);
    data([1:2,end-1:end],:) = [];
    data = reshape(data,datsiz);

fprintf('[%s] Averaged %.1f%% samples are filled with the most likely values\n', mfilename, mean(omitrat)*100);
else
fprintf('[%s] Averaged %.1f%% samples are filled with NaN\n', mfilename, mean(omitrat)*100);
end
