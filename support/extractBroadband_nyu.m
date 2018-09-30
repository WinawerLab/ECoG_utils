function broadband = extractBroadband_nyu(x, srate, method, bands)
%Compute time varying broadband envelope of a time series
% broadband = extractBroadband(x, srate, method, bands)
%
% Inputs
%   x:      data (time x n) (n is number of channels or epochs)
%
%   srate:  sample rate (Hz) [default = 1000]
%
%   method: method for computing broadband, can be a number 1-4, or a
%            string. [default = 1]
%               string
%
%                 1 'abs(hilbert(mean(whiten(bp(x)))))';
%                 2 'abs(hilbert(mean(whiten(bp(x))))).^2';
%                 3 'geomean(abs(hilbert(whiten(bp(x)))))';
%                 4 'geomean(abs(hilbert(whiten(bp(x))).^2)';
%                 5 'geomean(abs(hilbert(bp(x)).^2)';
%   bands:  required for methods 2-5. Can be a matrix (number of bands x 2)
%                     or a cell array {[lb, ub], width}
%                       [default = {[60 200], 20]}


if ~exist('srate', 'var')  || isempty(srate),  srate = 1000; end
if ~exist('method', 'var') || isempty(method), method = 1;   end
if ~exist('bands', 'var')  || isempty(bands),  bands = {[60 200], 20}; end


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
bp  = zeros([size(x) size(bands,1)], 'single');
for ii = 1:size(bands,1)
    bp(:,:, ii) = butterpass_eeglabdata_nyu(x,bands(ii,:),srate);
end

% if only one time series, then eliminate singleton dimension
if size(bp, 2) == 1, bp = squeeze(bp); end

% which dimension represents the multiple bands?
banddim = length(size(bp));

whiten = @(x) (x - mean(x(:)))./ diff(prctile(x, [.25 .75]));
%whiten = @(x) x;

switch method
    case {'abs(hilbert)', 1} 
        % -- Method 1: 'abs(hilbert(mean(whiten(bp(x)))))';  
        broadband = abs(hilbert(mean(whiten(bp),banddim)));
    
    case {'abs(hilbert(bp))', 2}
        % -- Method 2: 'abs(hilbert(mean(whiten(bp(x))))).^2';
        broadband = abs(hilbert(mean(whiten(bp), banddim))).^2;
    
    case {'abs(hilbert(sum(whiten(bp))))', 3}
        % -- Method 3: geomean(abs(hilbert(whiten(bp(x))))
        broadband = geomean(abs(hilbert(whiten(bp))),  banddim);

    case {'sum(abs(hilbert(whiten(bp))))', 4}
        % -- Method 4: geomean(abs(hilbert(whiten(bp(x)))).^2)
        broadband = geomean(abs(hilbert(whiten(bp))).^2, banddim); 
	
    case {'sum(abs(hilbert(bp)))', 5}
        % -- Method 5: geomean(abs(hilbert(bp(x))).^2)
        broadband = geomean(abs(hilbert(bp)).^2, banddim);  
end

return