function [coef,predictor] = ecog_regressERP(data_epoch,t_base,stims,stim4predict,maxlag,negregress)
% function to regress the ERP from some data
% USAGE:
% [coef,predictor] = ecog_regressERP(data_epoch,t_base,stims,stim4predict)
% [coef,predictor] = ecog_regressERP(data_epoch,t_base,stims,stim4predict,maxlag,negregress)
% 
% data_epoch % electrodes X time X epochs
% t_base % indices of the baseline period (ex: [1:100])
% stims % different code for different conditions
% stim4predict % the number of the condition code to construct the predictor
% maxlag % maximum time sample of lag for regression (default: 0)
% negregress % if false, adopt the lag when the coefficient is largest (default)
%              if true, adopt the lag when the absolute coefficient is larget
% 
% See also, ecog_regressout

% 20190903 Yuasa: modified from ecog_spectra.m in ECoG_utils
% 20200303 Yuasa: make computations effective & use 'parfor'
% 20210413 Yuasa: use same predictor for all stims
% 20210425 Yuasa: add maxlag option

%-- check data_epoch
narginchk(4,6);
if nargin < 6,  negregress = false;     end
if nargin < 5,  maxlag = 0;             end
dim_tim = size(data_epoch)==length(t_base);
dim_trl = size(data_epoch)==length(stims);
assert(any(dim_tim),'t_base is invalid');
assert(any(dim_trl),'stims is invalid');
if dim_tim(2), dim_tim = 2;
else,          dim_tim = find(dim_tim,1);
end
dim_trl(dim_tim) = false;
assert(any(dim_trl),'t_base or stims is invalid');
if dim_trl(3), dim_trl = 3;
else,          dim_trl = find(dim_trl,1);
end
dim_ch  = true(1,ndims(data_epoch));
dim_ch([dim_tim, dim_trl]) = false;
dim_ch = find(dim_ch,1);
assert(any(ismember(stim4predict,stims)),'stim4predict is invalid');

%-- permute data_epoch
data_epoch = permute(data_epoch,[dim_ch,dim_tim,dim_trl]);

%-- t_base must be logical array
if numel(t_base)==2 && diff(t_base)~=1
    warning('t_base might be invalid');
end

% %-- convert ordinal indices
% stim4predict = find(stims == stim4predict(1),1);
% stims = grp2idx(stims);
% stim4predict = stims(stim4predict);

%-- baseline correct
data_epoch = bsxfun(@minus, data_epoch, mean(data_epoch(:,t_base,:),2));

%-- regress erp out
predictor = zeros(size(data_epoch,[1,2]));
coef = nan(size(data_epoch,[1,3]));
for k=1:size(data_epoch,1)%channels
    disp(['regressing erp from channel ' int2str(k)])
    %-- make reference
    av_erp=zeros(size(data_epoch,2),1);
    av_erp(:,1) = mean(data_epoch(k,:,ismember(stims,stim4predict)),3,'omitnan');
    %-- regress ERP out
    if all(diff(av_erp)==0)     % skip if regressor is flat
        coef(k,:) = 0;
    elseif maxlag == 0          % regress without lag
        %%-- apply regress
        parfor m=1:size(data_epoch,3)%epochs
            x=reshape(data_epoch(k,:,m),[],1);
            if ~all(isnan(x))
                coef(k,m) = regress(x,av_erp);
            end
        end
    else                        % regress with lags
        parfor m=1:size(data_epoch,3)%epochs
            x=reshape(data_epoch(k,:,m),[],1);
            %%-- estimate time lag
            [cr,tlags] = xcorr(x,av_erp,maxlag);
            if all(diff(cr)==0)
                tlags = 0; dt = 1;   % avoid empty output for flat cr
            elseif negregress
                [~,dt] = max(abs(cr));
            else
                [~,dt] = max(cr);
            end
            %%-- apply regress
            reg_erp = circshift(av_erp,tlags(dt));
            coef(k,m) = regress(x,reg_erp);
        end
    end
    predictor(k,:,:) = av_erp;
end

%-- reverse permute data_epoch
rev_idx([dim_ch,dim_tim,dim_trl]) = 1:ndims(data_epoch);
predictor = permute(predictor,rev_idx);
