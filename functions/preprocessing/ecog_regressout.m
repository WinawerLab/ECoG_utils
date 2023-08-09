function [data_out,coefs,predictor,lags] = ecog_regressout(data_epoch,t_base,stims,maxlag,t_ref,negcoef4lag,stim4predict,usepar)

% function to regress out ERP from data
% USAGE:
% data_out = ECOG_REGRESSOUT(data_epoch,t_base,stims)
%   Regressor is computed for each stimulus based on 'stims'.
%   Baseline correction is applied before regression based on 't_base'.
% 
% data_out = ECOG_REGRESSOUT(data_epoch,t_base,stims,[maxlag,t_ref,negcoef4lag])
%   If 'maxlag' is specified, regression is computed allowing a lag in the range of [-maxlag, maxlag].
%   The best time lag for each epoch is estimated based on cross-correlation.
%   Time ranges to estimate coefficients of regression can limit in 't_ref'.
%   If 'negcoef4lag' is true, the lag is set to the best fit even the coefficient is negative. 
% 
% data_out = ECOG_REGRESSOUT(data_epoch,t_base,stims,maxlag,t_ref,negcoef4lag,stim4predict)
%   A common regressor is applied instead of stimulus-dependent regressors.
%   The common regressor is computed from the stimuli specified in 'stim4predict'.
% 
% [data_out,coef,predictor,lag] = ECOG_REGRESSOUT(data_epoch,t_base,stims [,maxlag,t_ref,negcoef4lag])
%   returns coefficients of regression, regressor, and estimated lag.
% 
% ECOG_REGRESSOUT(...,useparallel=false) suppress to use parallel computing.
% 
% Inputs: 
%   data_epoch   % electrodes X time X epochs
%   t_base       % logical array to indicate the baseline time points
%   stims        % different code for different conditions
%   maxlag       % allows lag in the range from -maxlag to maxlag for regression
%                  if 0, regression is applied without lag (default: [0])
%                  (ex: [5] for maximum 10 ms lag at 500 Hz sampling rate)
%   t_ref        % logical array to indicate the reference time points to apply regression (default: all time points)
%   negcoef4lag  % if false, adopt the lag when the coefficient is largest (default)
%                  if true, adopt the lag when the absolute coefficient is largest
%   stim4predict % array of stimulus number indicated in 'stims'
%                  if not empty, a common regressor averaging 'stim4predict' epochs is used instead
% 
% Outputs: 
%   data_out     % data_epoch after regressing out the EPR
%   coef         % coefficients of regressor
%   predictor    % regressor (ERP)
%   lag          % lag where regression is applied
% 
% See also, ecog_regressERP, regress, xcorr

% 20190903 Yuasa: modified from ecog_spectra.m in ECoG_utils
% 20200303 Yuasa: make computations effective & use 'parfor'
% 20210810 Yuasa: merge ecog_regresslag into this function
% 20220420 Yuasa: enable to output coef and predictor
% 20220825 Yuasa: add stim4predict option
% 20230720 Yuasa: add useparallel option

narginchk(3,inf);
%-- check data_epoch
dim_tim = size(data_epoch)==length(t_base);
dim_trl = size(data_epoch)==length(stims);
assert(any(dim_tim),'t_base is invalid');
assert(any(dim_trl),'stims is invalid');
if dim_tim(2), dim_tim = 2;
else,          dim_tim = find(dim_tim,1);
end
dim_trl(dim_tim) = false;
assert(any(dim_trl),'t_base or stims is invalid');
if length(dim_trl)>=3 && dim_trl(3), dim_trl = 3;
else,          dim_trl = find(dim_trl,1);
end
dim_ch  = true(1,max(3,ndims(data_epoch)));
dim_ch([dim_tim, dim_trl]) = false;
dim_ch = find(dim_ch,1);

%-- permute data_epoch
data_epoch = permute(data_epoch,[dim_ch,dim_tim,dim_trl]);

%-- get inputs
nsample  = size(data_epoch,2);
if nargin < 8
    usepar = true;
end
if nargin < 7
    stim4predict = [];
end
if nargin < 6 || isempty(negcoef4lag)
    negcoef4lag = false;
end
if nargin < 5 || isempty(t_ref)
    t_ref = true(1,nsample);
end
if nargin < 4 || isempty(maxlag)
    maxlag = 0;
end

%-- t_base must be logical array
if numel(t_base)==2 && diff(t_base)>1
    warning('t_base might be invalid');
end
if numel(t_ref)==2 && diff(t_ref)>1
    warning('t_ref might be invalid');
end

%-- convert ordinal indices
[stims,~,grplst] = grp2idx(stims);
if ~isempty(stim4predict)    % common regressor
    %-- check stim4predict
    stim4predict = find(ismember(grplst,stim4predict));
    assert(~isempty(stim4predict),'stim4predict is invalid');
end

%-- baseline correct
if any(t_base)
    data_epoch = bsxfun(@minus, data_epoch, mean(data_epoch(:,t_base,:),2,'omitnan'));
end

%-- start parallel pool
poolsize = 0;
if usepar && ~isempty(gcp)
  poolsize = gcp;  poolsize = poolsize.NumWorkers;
end

%-- regress erp out
data_out  = nan(size(data_epoch));
predictor = zeros(size(data_epoch));
coefs     = nan(size(data_epoch,[1,3]));
lags      = zeros(size(data_epoch,[1,3]));
for k=1:size(data_epoch,1)%channels
    disp(['regressing erp from channel ' int2str(k)])
    %-- make reference
    if isempty(stim4predict)    % basic 
        av_erp_list = zeros(nsample,max(stims(:)));
        for m = unique(stims(:)')
            av_erp_list(:,m) = mean(data_epoch(k,:,stims==m),3,'omitnan');
        end
    else                        % common regressor
        av_erp_list = mean(data_epoch(k,:,ismember(stims,stim4predict)),3,'omitnan');
        av_erp_list = repmat(av_erp_list(:),1,max(stims(:)));
    end
    %-- regress ERP out
    parfor (m=1:size(data_epoch,3),poolsize) %epochs
      av_erp            = av_erp_list(:,stims(m));
      av_erp_pad        = zeros(nsample,1);
      av_erp_pad(t_ref) = av_erp(t_ref,:);
      x=reshape(data_epoch(k,:,m),nsample,1);
      if all(diff(av_erp_pad)==0)    % skip if regressor is flat
        data_out(k,:,m)=x;
        coef = 0;
      elseif maxlag == 0             % regress without lag
        %%-- apply regress
        if all(isnan(x))
            coef = 0;
        else
            coef = regress(x,av_erp_pad);
        end
        data_out(k,:,m) = x - coef .* av_erp;
      else                          % regress with lag
        x_pad        = zeros(nsample,1);
        x_pad(t_ref) = x(t_ref);
        % estimate time lag
        [cr,tlags] = xcorr(x_pad,av_erp_pad,maxlag);
        if all(diff(cr)==0)
          tlags = 0; dt = 1;   % avoid empty output for flat cr
        elseif negcoef4lag
            [~,dt] = max(abs(cr));
        else
            [~,dt] = max(cr);
        end
        % regress ERP out
        reg_erp     = circshift(av_erp,tlags(dt));
        reg_erp_pad = circshift(av_erp_pad,tlags(dt));
        coef = regress(x,reg_erp_pad);
        data_out(k,:,m) = x - coef .* reg_erp;
        lags(k,m)  = tlags(dt);
      end
      coefs(k,m) = coef;
    end
    predictor(k,:,:) = av_erp_list(:,stims);
end

%-- reverse permute data_epoch
rev_idx([dim_ch,dim_tim,dim_trl]) = 1:ndims(data_out);
data_out  = permute(data_out,rev_idx);
if ~isempty(stim4predict)    % common regressor
    predictor = predictor(:,:,1);
end
predictor = permute(predictor,rev_idx);
