function [data,omits] = ecog_omitNoisySamples(data,srate,varargin)

% data = ecog_omitNoisySamples(data, srate, noise_freq, [nharmonic, alpha, noise_mindur, noise_expand, noise_ratio])
%   fills NaN values in time points where artifact noise is significantly large.
%
% [data, omits] = ecog_omitNoisySamples(data, srate, noise_freq, ...)
%   returns 'omits', a logical array indicating omitted samples.
%
% [data, omits] = ecog_omitNoisySamples(..., 'simple')
%   detects outliers only based on the amplitude.
%
% data = ecog_omitNoisySamples(..., 'fill')
%   fills missing values after omitting noisy samples.
%
% INPUTS:
%   data:          (time X channels)
%   srate:         Sampling rate (Hz)
%   noise_freq:    Artifact noise frequency (Hz) [default = 60]
%   nharmonic:     Nth harmonics to be notched [default = 3]
%   alpha:         Exclusion threshold for noise detection (scalar) [default = 0.003 (0.3%)]. 
%                  Can be a vector, [threshold, duration of segments (s)],
%                  noise detection is performed in segments with at least 20% overlap,
%                  and results are merged conservatively.
%   signal_mindur: Signal between noise must appear for at least S seconds [default = 0.50 <500ms>]
%   noise_mindur:  Noise must appear for at least S seconds [default = 0.02 <20ms>]
%   noise_ratio:   X*100% samples in the noise_mindur need to be detected as noise [default = 0.8 <80%>]
%   noise_expand:  Expand noise detected samples for before and after S seconds [default = 0.10 <100ms>]

% 20230816 Yuasa
% 20231105 Yuasa: embeded segmentation mode

%-- Set parameters
narginchk(2,inf);
%--- extra options
extropts = cellfun(@(v) ischar(v)||isstring(v),varargin);
issimple = ismember('simple',varargin(extropts));
isfill   = ismember('fill',varargin(extropts));
%--- basic options
[noise_freq, nharmonic,alpha,signal_mindur,noise_mindur,noise_ratio,noise_expand] = parseInputs(varargin{~extropts});
if isempty(nharmonic), nharmonic = 3; end      % Nth-harmonic
if isempty(alpha) || isnan(alpha(1)), alpha(1) = 0.003; end    % [0,1]
if isempty(signal_mindur), signal_mindur = 0.50; end   % [seconds]
if isempty(noise_mindur), noise_mindur = 0.02; end   % [seconds]
if isempty(noise_ratio), noise_ratio = 0.8; end   % [0,1]
if isempty(noise_expand), noise_expand = 0.10; end   % [seconds]
%--- set up options
nharmonic      = round(nharmonic(1));
signal_minsmpl = max(round(signal_mindur(1)*srate),1);
noise_minsmpl  = max(round(noise_mindur(1)*srate),1);
noise_ratio    = noise_ratio(1);
noise_expsmpl  = max(round(noise_expand(1)*srate),1);
datleng        = size(data,1);
if isscalar(alpha)
    segment_smpl = datleng;
else
    segment_smpl = min(datleng,round(alpha(2)*srate));
    alpha        = alpha(1);
end
assert(alpha>0&alpha<1, 'Invalid alpha is specified.');

if issimple
%-- Simple setting
segment_st = 1;
else
%-- Segments
segment_minovl = 0.2;   % 20% overlap
segment_n = max( 1, (datleng-segment_smpl.*segment_minovl) ./ (segment_smpl.*(1-segment_minovl)) );
if (segment_n - fix(segment_n)) < 0.05
    segment_n = fix(segment_n);
    segment_smpl = ceil( datleng ./ (segment_n.*(1-segment_minovl) + segment_minovl) );
else
    segment_n = ceil(segment_n);
end
segment_st = round(linspace(1, datleng-segment_smpl+1, segment_n));

%-- Expand notch frequencies to include all harmonics
if ~isempty(noise_freq)
frange = 5;
noise_freq = reshape(noise_freq(:).*(1:nharmonic),[],1) + [-1 1].*frange;
noise_freq(noise_freq>srate/2) = srate/2;
noise_freq(diff(noise_freq,1,2)==0,:) = [];
if size(noise_freq,1)>1
refer_freq = round(10.^movmean(log10(noise_freq),2,1,EndPoint='discard'));
else
refer_freq = noise_freq + 10*frange;
end
end

%-- Compute power timecourse at noise frequnecy
fprintf('[%s] Identifying and omitting noisy samples...\n', mfilename);
if ~isempty(noise_freq)
[noise_signal] = ecog_extractBroadband(data, srate, [], noise_freq);
[refer_signal] = ecog_extractBroadband(data, srate, [], refer_freq);
end
[std_signal]   = movstd(data,max(noise_minsmpl*3,100),0,1);
end

%-- Fill NaN for each channel
omits   = false(size(data));
for ich = 1:prod(size(data,2:ndims(data)))
    outl_smpls  = false(datleng,1);
    if issimple
    noise_smpls = false(datleng,1);
    else
    noise_smpls = true(datleng,1);
    for segstidx = segment_st
        segidx = (1:segment_smpl)-1+segstidx;
        %-- Detect sharp noise samples
        trgt_samples = abs(diff([mean(data(segidx,ich),'omitnan');data(segidx,ich)]));
        outl_smplseg   = detectNoise(trgt_samples,alpha,'Gamma',false);
        %-- Detect flat signal
        trgt_samples = 1./std_signal(segidx,ich).^2;  trgt_samples(outl_smplseg) = nan;
        noise_smplseg = ...
                        detectNoise(trgt_samples,alpha,'Loglogistic',true,noise_minsmpl,noise_ratio);
        if ~isempty(noise_freq)
        %-- Detect absolutely significant large noise
        trgt_samples = sqrt(noise_signal(segidx,ich));  trgt_samples(outl_smplseg) = nan;
        noise_smplseg = noise_smplseg | ...
                        detectNoise(trgt_samples,alpha,'Loglogistic',true,noise_minsmpl,noise_ratio);
        %-- Detect relatively significant large noise
        trgt_samples = sqrt(noise_signal(segidx,ich)./refer_signal(segidx,ich));  trgt_samples(outl_smplseg) = nan;
        noise_smplseg = noise_smplseg | ...
                        detectNoise(trgt_samples,alpha,'Loglogistic',true,noise_minsmpl,noise_ratio);
        end
        outl_smpls(segidx)   = outl_smpls(segidx) | noise_smplseg;
        noise_smpls(segidx) = noise_smpls(segidx) & noise_smplseg;
    end
    end
    %-- Detect outlier
    trgt_samples = abs(data(:,ich)); trgt_samples(outl_smpls|noise_smpls) = nan;
    outl_smpls = outl_smpls | ...
                 trgt_samples > min(prctile(trgt_samples,[90 99 99.9 99.99]).*[100 10 5 2]);
    %-- Apply threshold on noise_samples & Omit short signal
    noise_smpls = outl_smpls | ...
        threshNoise(outl_smpls|noise_smpls,0,noise_ratio,noise_expsmpl,signal_minsmpl);
    %-- Fill NaN
    data(noise_smpls,ich) = nan;
    omits(:,ich)            = noise_smpls;
end
omitrat = sum(omits,1)./size(omits,1);
fprintf('[%s] Averaged %.1f%% of samples are identified as noise\n', mfilename, mean(omitrat)*100);

%-- Filling missing
if isfill
    datsiz = size(data);
    padding = zeros([2,size(data,2:ndims(data))]);
    data = cat(1,padding,data,padding);
    data = fillmissing(data,'pchip',1);
    data([1:2,end-1:end],:) = [];
    data = reshape(data,datsiz);

    fprintf('[%s] Noisy samples are filled with the most likely values\n', mfilename);
end

end


%%% Sub function
function varargout = parseInputs(varargin)
%-- passthrough inputs 
varargout = cell(1,nargout);
varargout(1:min(nargout,nargin)) = varargin(1:min(nargout,nargin));

end


function noise_samples = detectNoise(noise_signal,alpha,distname,updatealpha,noise_minsmpl,noise_ratio,noise_expsmpl)

if ~exist('distname','var') || isempty(distname), distname = 'Loglogistic';  end
if ~exist('updatealpha','var') || isempty(updatealpha), updatealpha = false;  end
if ~exist('noise_minsmpl','var'), noise_minsmpl = [];  end
if ~exist('noise_ratio','var'),   noise_ratio = [];    end
if ~exist('noise_expsmpl','var'), noise_expsmpl = [];  end
%-- Statistical distoribution fitting
if min(noise_signal,[],'all','omitnan')>=0
noise_signal = noise_signal ./ median(noise_signal,'omitnan');
end
estdist = fitdist(noise_signal,distname);
%-- Detect noise
noise_samples = noise_signal > estdist.icdf(1-alpha);
%-- Update alpha if too much samples are detected (over 3*alpha%)
alpha_orig = alpha; iratio = 3;
while updatealpha && sum(noise_samples)/sum(~isnan(noise_signal)) > alpha_orig*iratio
    alpha = alpha/2; iratio = iratio+0.5;
    noise_samples = noise_signal > estdist.icdf(1-alpha);
end

%-- Apply threshold
if ~isempty(noise_minsmpl)|| ~isempty(noise_expsmpl)
noise_samples = threshNoise(noise_samples,noise_minsmpl,noise_ratio,noise_expsmpl);
end

end


function noise_samples = threshNoise(noise_samples,noise_minsmpl,noise_ratio,noise_expsmpl,signal_minsmpl)
if ~exist('noise_minsmpl','var')  || isempty(noise_minsmpl)  || noise_minsmpl<1,  noise_minsmpl = 1;   end
if ~exist('noise_ratio','var')    || isempty(noise_ratio)    || noise_ratio>1,    noise_ratio = 1;     end
if ~exist('noise_expsmpl','var')  || isempty(noise_expsmpl)  || noise_expsmpl<1,  noise_expsmpl = 1;   end
if ~exist('signal_minsmpl','var') || isempty(signal_minsmpl) || signal_minsmpl<1, signal_minsmpl = 1;  end
%-- Apply minimum length & max ratio requirement
noise_samples = movmax(conv(noise_samples,ones(noise_minsmpl,1)/noise_minsmpl)...
                            ,noise_minsmpl,EndPoints='discard')...
                        >= noise_ratio;
%-- Expand noise samples
noise_samples = movmean(noise_samples,2*noise_expsmpl) > 0;
%-- Omit short signal
noise_samples = movmin(conv(noise_samples,ones(signal_minsmpl,1)/signal_minsmpl)...
                            ,signal_minsmpl,EndPoints='discard')...
                        > 0;

% %%% precise onset & offset
% noise_gapsmpl = round(noise_minsmpl .* (1-noise_ratio));
% %--- Ignore small gap 
% noise_clstr   = cell2mat(arrayfun(@(t)t.*ones(t,1),diff(find([1;diff(noise_samples);1])),'UniformOutput',false));
% noise_samples(noise_clstr.*~noise_samples <= noise_gapsmpl) = true;
% %--- Ignore short noise 
% noise_clstr   = cell2mat(arrayfun(@(t)t.*ones(t,1),diff(find([1;diff(noise_samples);1])),'UniformOutput',false));
% noise_samples(noise_clstr.*noise_samples < noise_minsmpl) = false;

end
