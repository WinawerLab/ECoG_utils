function [band_sig,norder]=butterpass_zerophase(signal,band,srate,norder,decreaseorder)
%% Bandpass or lowpass filter a time series using a zero-phase butterworth filter
% [band_sig,order]=butterpass_zerophase(signal,band,srate,order)
%
% Written in the same rule as butterpass_eeglabdata_nyu
%
% INPUTS:
%   signal: ( time X channels )    
%   band:   Cutoff frequency (Hz)
%           [Fhp Flp] for bandpass filter
%           [Flp] for lowpass filter
%   srate:  sampling rate (Hz) [default = 1000]
%   order:  Nth-order filter [default = 4]
%
% OUTPUTS: 
%   band_sig ( time X channels)  : band-passed signal
%
% Example 1: visualize a filter
%   butterpass_zerophase([],[1 60]);
% Example 2: lowpass filter a noisy vector
%  x = randn(1000,1);
%  bp = butterpass_zerophase(x,[60], 1000, 4);
%  plot(1:length(x), x, 1:length(bp), bp)
% Example 3: bandpass a noisy matrix
%  x = randn(1000,10);
%  bp = butterpass_zerophase(x,[10 60], 1000, 4);
%  plot(1:length(x), x, 'r', 1:length(bp), bp, 'k')
% 
% See also FILTFILT, BUTTERPASS_EEGLABDATA_NYU

% 20220801 Yuasa

% If there is no signal, then plot the filter and return
if ~exist('signal', 'var') || isempty(signal) 
    plotFilter = true; 
else, plotFilter = false;
end

if ~exist('srate', 'var') || isempty(srate), srate = 1000; end   % sample rate
if ~exist('norder', 'var') || isempty(norder), norder = 4; end   % Nth-order
if ~exist('decreaseorder', 'var') || isempty(decreaseorder), decreaseorder = true; end   % allow to modify order if output is invalid


% Design Butterworth IIR filter.
if isscalar(band)
  Hd = designfilt('lowpassiir', 'FilterOrder', norder, ...
       'HalfPowerFrequency', band,...
       'SampleRate', srate,'DesignMethod','butter');
else
  Hd = designfilt('bandpassiir','FilterOrder',norder*2, ...
         'HalfPowerFrequency1',min(band),'HalfPowerFrequency2',max(band), ...
         'SampleRate',srate,'DesignMethod','butter');
end

% visualize
if plotFilter
    fvtool(Hd)
    band_sig = [];
    return
end

% Apply filter
band_sig = filtfilt(Hd, signal); 

% Check output
normpow = rms(signal,1);
normpow = mean(normpow - rms(band_sig,1),'all','omitnan') ./  mean(normpow,'all','omitnan');
if normpow < -0.2
    if norder > 1 && decreaseorder
        warning('The size of filtering oder is not appropriate. Trying order=%d.',norder-1);
        [band_sig,norder]=butterpass_zerophase(signal,band,srate,norder-1);
    else
        warning('The size of filtering oder is not appropriate.');
    end
end









