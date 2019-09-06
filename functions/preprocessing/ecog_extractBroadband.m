function [broadband, methodstr, bands] = ecog_extractBroadband(x, srate, method, bands)
% Compute time varying broadband envelope of a time series
% broadband = extractBroadband(x, srate, method, bands)
%
% Inputs
%   x:      data (time x n) (n is number of channels or epochs)
%
%   srate:  sample rate (Hz) [default = 1000]
%
%   method: a function handle that specifies the method for computing 
%           broadband that takes as input band-passfiltered data (bp) 
%           of dimensions (number of bands x time x n). 
%           [default: @(bp) geomean(abs(hilbert(bp)).^2)]
%            
%   bands:  a matrix (number of bands x 2) or a cell array {[lb, ub], width}
%           [default = {[70 200], 20}]
% 
% Note: the bandpass-filtered data is stacked in the first dimension of bp,
% and the mean (if specified in the function handle) is taken over that
% dimension. If filtering a single band, the method function handle should
% not include an averaging step, or should explicitly specify that the mean
% is to be taken over the first, singleton dimension; see Example 3 below.
%
% Example 1: default settings:
%   data = randn(10000,1);  
%   [broadband, methodstr] = ecog_extractBroadband(data);
%   figure, 
%   subplot(2,1,1); plot(1:length(data), data)
%   subplot(2,1,2); plot(1:length(data), broadband);
%
% Example 2: compare different bandwidths:
%   data = randn(10000,1);  
%   srate = 1000;
%   method = @(bp) geomean(abs(hilbert(bp)).^2);
%   bands  = {[60 200], 35};
%   [broadband1, methodstr1] = ecog_extractBroadband(data);
%   [broadband2, methodstr2] = ecog_extractBroadband(data, [], [], bands);
%   figure, 
%   subplot(2,1,1); plot(1:length(data), data)
%   subplot(2,1,2); plot(1:length(data), broadband1, 1:length(data), broadband2);
% 
%   Example 3: use a single band and adapt method accordingly:
%   data = randn(10000,1);  
%   bands  = {[60 200], 140};
%   method = @(bp) abs(hilbert(bp)).^2;
%   [broadband, methodstr] = ecog_extractBroadband(data, [], method, bands);
%   figure, 
%   subplot(2,1,1); plot(1:length(data), data)
%   subplot(2,1,2); plot(1:length(data), broadband);

if ~exist('srate', 'var')  || isempty(srate),  srate = 1000; end
if ~exist('bands', 'var')  || isempty(bands),  bands = {[60 200], 20}; end
if ~exist('method', 'var') || isempty(method), method = @(bp)geomean(abs(hilbert(bp)).^2); end

% format the bands input
if isa(bands, 'cell')
    
    % Entire range for broadband
    band_rg  = bands{1}; 
    
    % Bin width 
    band_w   = bands{2}; 

    % All bins
    lb       = band_rg(1):band_w:band_rg(2)-band_w;
    ub       = lb+band_w;
    bands   = [lb; ub]';
end
    
% band pass filter each sub-band
% the first dimension represents the multiple bands
bp  = zeros([size(bands,1) size(x)], 'double');

for ii = 1:size(bands,1)
    fprintf('[%s] Filtering signal in band %d-%d\n',mfilename,bands(ii,1),bands(ii,2));
    bp(ii,:,:) = butterpass_eeglabdata_nyu(x,bands(ii,:),srate); 
    % divide power by band width?
    %tmp = butterpass_eeglabdata_nyu(x,bands(ii,:),srate);
    %tmp = tmp / diff(bands(ii,:)); 
    %bp(ii,:,:) = tmp;
end

% if only one time series, then eliminate singleton dimension
% if size(bp, 2) == 1, bp = squeeze(bp); end

%banddim = length(size(bp));

%whiten = @(x) (x - mean(x(:)))./ diff(prctile(x, [.25 .75]));
%whiten = @(x) x;

fprintf('[%s] Computing broadband... \n', mfilename);

broadband = method(bp);
methodstr = func2str(method);

% if only one band is used, eliminate singleton dimension
if size(broadband, 1) == 1, broadband = squeeze(broadband); end

% old method specification:
%
% switch method
%     case {1}
%         % -- Method 1: mean before hilbert; regular mean; whiten; amplitude
%         broadband = abs(hilbert(mean(whiten(bp),banddim)));
%         methodstr = 'abs(hilbert(mean(whiten(bp))))';
%         
%     case {2}
%         % -- Method 2: mean before hilbert; regular mean; whiten; power
%         broadband = abs(hilbert(mean(whiten(bp), banddim))).^2;
%         methodstr = 'abs(hilbert(mean(whiten(bp)))).^2';
%         
%     case {3}
%         % -- Method 3: mean after hilbert; geomean; whiten; amplitude
%         broadband = geomean(abs(hilbert(whiten(bp))),  banddim);
%         methodstr = 'geomean(abs(hilbert(whiten(bp)))';
%         
%     case {4}
%         % -- Method 4: mean after hilbert; geomean; whiten; power
%         broadband = geomean(abs(hilbert(whiten(bp))).^2, banddim);
%         methodstr = 'geomean(abs(hilbert(whiten(bp))).^2)';
%     
%     case {5}
%         % -- Method 5: mean after hilbert; geomean; amplitude
%         broadband = geomean(abs(hilbert(bp)), banddim);
%         methodstr = 'geomean(abs(hilbert(bp)))';
%         
%     case {6}
%         % -- Method 6: mean after hilbert; geomean; power
%         %--> BEST IN SIMULATION
%         broadband = geomean(abs(hilbert(bp)).^2, banddim);
%         methodstr = 'geomean(abs(hilbert(bp)).^2)';
%     
%     case {7}
%         % -- Method 7: mean after hilbert; regular mean; whiten; power
%         broadband = mean(abs(hilbert(whiten(bp))).^2, banddim);
%         methodstr = 'mean(abs(hilbert(whiten(bp))).^2)';
%        
%     case {8}
%         % -- Method 8: mean after hilbert; regular mean; logpower
%         broadband = mean(log10(abs(hilbert(bp)).^2), banddim);
%         methodstr = 'mean(log10(abs(hilbert(bp)).^2))';
%     otherwise
%         error('Unknown broadband method %s', method);
% end
% broadband = broadband';

%fprintf('done\n');

return

