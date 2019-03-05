function [band_sig]=butterpass_eeglabdata_nyu(signal,band,srate, Rp, Rs, bw)
%% Bandpass filter a time series using a butterworth filter
% [band_sig]=butterpass_eeglabdata(signal,band,srate, Rp, Rs, bw)
%
% Modified from Dora's code by Jon
%
% INPUTS:
%   signal: ( time X channels)    
%   band:   [start:stop] or [start stop]
%   srate:  sampling rate (Hz) [default = 1000]
%   Rp:     Passband ripple in decibels [default = 3]
%   Rs:     Stopband attenuation in decibels [default = 60]
%   bw:     Stop band is band(1)*bw and band(2)/bw [default = .5]
%
% OUTPUTS: 
%   band_sig ( time X channels)  : band-passed signal
%
% Example 1: visualize a filter
%   butterpass_eeglabdata([],[80 120]);
% Example 2: filter a noisy vector
%  x = randn(1000,1);
%  bp = butterpass_eeglabdata(x,[80 100], 1000, 3, 60, .5);
%  plot(1:1000, x, 1:1000, bp)
% Example 3: filter a noisy matrix
%  x = randn(1000,10);
%  bp = butterpass_eeglabdata(x,[80 100], 1000, 3, 60, .5);
%  plot(1:1000, x, 'r', 1:1000, bp, 'k')

% If there is no signal, then plot the filter and return
if ~exist('signal', 'var') || isempty(signal) 
    plotFilter = true; 
else, plotFilter = false;
end

if ~exist('srate', 'var') || isempty(srate), srate = 1000; end   % sample rate
if ~exist('Rp', 'var') || isempty(Rp), Rp = 3; end   % Passband ripple in decibels
if ~exist('Rs', 'var') || isempty(Rs), Rs = 60; end  %  Stopband attenuation in decibels
if ~exist('bw', 'var') || isempty(bw), bw = .5; end   % Stop band is band(1)*bw and band(2)/bw 


% Butterworth filter that loses no more than Rp dB ripple in the passband
% and has at least Rs dB of attenuation in the stopband.

Fstop1 = band(1)*bw;  % First Stopband Frequency
Fpass1 = band(1);     % First Passband Frequency
Fpass2 = band(end);   % Second Passband Frequency
Fstop2 = band(end)/bw;% Second Stopband Frequency
Astop1 = Rs;          % First Stopband Attenuation (dB)
Apass  = Rp;          % Passband Ripple (dB)
Astop2 = Rs;          % Second Stopband Attenuation (dB)
match  = 'stopband';  % Band to match exactly

nyqLimit = .5 * srate;
Fstop2 = min(nyqLimit, Fstop2);

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, ...
                      Astop2, srate);
Hd = design(h, 'butter', 'MatchExactly', match);

% visualize
if plotFilter
    fvtool(Hd, 'Fs', srate)
    band_sig = [];
    return
end


% measure time shift of filter
[gd,f] = grpdelay(Hd, 512 ,srate);
shift_frames = round(mean(gd(f>=band(1) & f <= band(2))));

% Band pass
band_sig = filter(Hd, signal); 

% correct for time shift of filter
band_sig = circshift(band_sig, -shift_frames,1);




