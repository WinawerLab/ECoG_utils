function [broadband, methodstr, bands] = ecog_extractBroadband(x, srate, method, bandopts)
% Compute time varying broadband envelope of a time series
% broadband = ecog_extractBroadband(x, srate, method, bandsopts)
%
% Inputs
%   x:      data (time x n) (n is number of channels or epochs)
%
%   srate:  sample rate (Hz) [default = 1000]
%
%   method: a function handle that specifies the method for computing 
%           broadband that takes as input bandpass filtered data (bp)  
%           of dimensions (time x n x number of bands). 
%           [default: @(bp,banddim) geomean(abs(hilbert(bp)).^2,banddim)]
%            
%   bandopts: a matrix (number of bands x 2) or a cell array {[lb, ub], width}
%           [default = {[70 200], 20}]
% 
% Note: If filtering a single band, the averaging step and therefore the
% banddim input argument can be eliminated (see Example 4 below).
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
%   bandopts  = {[60 200], 35};
%   [broadband1] = ecog_extractBroadband(data);
%   [broadband2] = ecog_extractBroadband(data, [], [], bandopts);
%   figure, 
%   subplot(2,1,1); plot(1:length(data), data)
%   subplot(2,1,2); plot(1:length(data), broadband1, 1:length(data), broadband2);
%
% Example 3: compare different methods:
%   data = randn(10000,1);  
%   srate = 1000;
%   method = @(bp,banddim) geomean(abs(hilbert(bp)),banddim); % take amplitude instead of power
%   [broadband1, methodstr1] = ecog_extractBroadband(data);
%   [broadband2, methodstr2] = ecog_extractBroadband(data, [], method);
%   figure, 
%   subplot(2,1,1); plot(1:length(data), data); 
%   subplot(2,1,2); plot(1:length(data), broadband1, 1:length(data), broadband2);
%   legend(methodstr1, methodstr2);
% 
% Example 4: use a single band and adapt method accordingly:
%   data = randn(10000,1);  
%   bandopts  = {[60 200], 140};
%   method = @(bp) abs(hilbert(bp)).^2;
%   [broadband, methodstr] = ecog_extractBroadband(data, [], method, bandopts);
%   figure, 
%   subplot(2,1,1); plot(1:length(data), data)
%   subplot(2,1,2); plot(1:length(data), broadband);

if ~exist('srate', 'var')    || isempty(srate),     srate = 1000; end
if ~exist('bandopts', 'var') || isempty(bandopts),  bandopts = {[60 200], 20}; end
if ~exist('method', 'var')   || isempty(method),    method = @(bp,banddim) geomean(abs(hilbert(bp)).^2,banddim); end

% format the bands input
if isa(bandopts, 'cell')
    
    % Entire range for broadband
    band_rg  = bandopts{1}; e
    
    % Bin width 
    band_w   = bandopts{2}; 

    % All bins
    lb       = band_rg(1):band_w:band_rg(2)-band_w;
    ub       = lb+band_w;
    bands   = [lb; ub]';
else
    bands = bandopts;
end
    
% band pass filter each sub-band
bp  = repmat({zeros(size(x), 'double')},size(bands,1),1);

% replace any nans in the data with zeros
if any(isnan(x(:)))
    warning('[%s] Found nans in the data! Replacing with zeros. Please check data for extensive nans \n', mfilename);
    x(isnan(x)) = 0;
end

% bandpass filter 
for ii = 1:size(bands,1)
    fprintf('[%s] Filtering signal in band %d-%d\n', mfilename, bands(ii,1),bands(ii,2));
    bp{ii} = butterpass_eeglabdata_nyu(x,bands(ii,:),srate); 
    % divide power by band width?
    %tmp = butterpass_eeglabdata_nyu(x,bands(ii,:),srate);
    %tmp = tmp / diff(bands(ii,:)); 
    %bp(ii,:,:) = tmp;
end
bp = cat(ndims(x)+1,bp{:});

% if only one time series, then eliminate singleton dimension
if size(bp, 2) == 1, bp = squeeze(bp); end

% compute broadband 
fprintf('[%s] Computing broadband... \n', mfilename);

methodstr = func2str(method);

if ~contains(methodstr, 'banddim')
    broadband = method(bp);
else
    banddim = ndims(bp);
    broadband = method(bp,banddim);
end

% old method specification:
%
%whiten = @(x) (x - mean(x(:)))./ diff(prctile(x, [.25 .75]));
%whiten = @(x) x;
%
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

