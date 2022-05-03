function [data_out,coefs,predictor,lags] = ecog_regressout(data_epoch,t_base,stims,maxlag,t_ref,negregress)

% function to regress out the ERP from some data
% USAGE:
% data_out = ECOG_REGRESSOUT(data_epoch,t_base,stims)
% data_out = ECOG_REGRESSOUT(data_epoch,t_base,stims [,maxlag,t_ref,negregress])
%   with specifying 'maxlag', it estimates time lag for regression based on cross-correlation
% 
% [data_out,coef,predictor,lag] = ECOG_REGRESSOUT(data_epoch,t_base,stims [,maxlag,t_ref,negregress])
% 
% Inputs: 
%   data_epoch % electrodes X time X epochs
%   t_base     % indices of the baseline period (ex: [1:100])
%   stims      % different code for different conditions
%   maxlag     % maximum time sample of lag for regression (default: 0)
%   t_ref      % indices of the reference period to apply regression (ex: [1:30])(default: all time points)
%   negregress % if allow negative regression during the time lag evaluation (default: false)
% 
% Outputs: 
%   data_out   % data_epoch after regressing out the EPR
%   coef       % coefficients of regressor
%   predictor  % regressor (ERP)
%   lag        % regressor (ERP)
% 
% See also, ecog_regressERP

% 20190903 Yuasa: modified from ecog_spectra.m in ECoG_utils
% 20200303 Yuasa: make computations effective & use 'parfor'
% 20210810 Yuasa: merge ecog_regresslag into this function
% 20220420 Yuasa: enable to output coef and predictor

narginchk(3,6);
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
if nargin < 6 || isempty(negregress)
    negregress = false;
end
if nargin < 5 || isempty(t_ref)
    t_ref = 1:nsample;
end
if nargin < 4 || isempty(maxlag)
    maxlag = 0;
end

%-- t_base must be logical array
if numel(t_base)==2 && diff(t_base)~=1
    warning('t_base might be invalid');
end
if numel(t_ref)==2 && diff(t_ref)~=1
    warning('t_ref might be invalid');
end

%-- convert ordinal indices
stims = grp2idx(stims);

%-- baseline correct
data_epoch = bsxfun(@minus, data_epoch, mean(data_epoch(:,t_base,:),2));

%-- regress erp out
data_out  = nan(size(data_epoch));
predictor = zeros(size(data_epoch));
coefs     = nan(size(data_epoch,[1,3]));
lags      = zeros(size(data_epoch,[1,3]));
for k=1:size(data_epoch,1)%channels
    disp(['regressing erp from channel ' int2str(k)])
    %-- make reference
    av_erp_list=zeros(size(data_epoch,2),max(stims(:)));
    for m=unique(stims(:)')
        av_erp_list(:,m) = mean(data_epoch(k,:,stims==m),3,'omitnan');
    end
    %-- regress ERP out
    parfor m=1:size(data_epoch,3)%epochs
      av_erp            = av_erp_list(:,stims(m))
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
        elseif negregress
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
predictor = permute(predictor,rev_idx);
