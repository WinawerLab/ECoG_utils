function [f,data_epoch_spectra,data_epoch] = ...
    ecog_crossspectra(data_epoch,stims,fft_w,fft_t,fft_ov,srate,reg_erp,f_lim,outtype,calcmode)

% Function to calculate cross-powerspectrum for across several ECoG
% channels and epochs using Matlab's cpsd function.
% 
% Input:
% data_epoch: channels X epochs X time
%
% stims: if regress_erp==1, the ERP calculated per condition can be regressed out.
%   The condition is indicated in stims, a vector with one value for each
%   epoch. If choosing to regress erp out, make sure data are baseline
%   corrected first, so the ERP can be calculated correctly!
% 
% fft_w: window size to devide the data_epoch into segment
% fft_t: time point to apply frequency analysis
% fft_ov: overlap size in each segment
% f_limits: output sepctra in [f_limits(1) f_limits(2)]
%
% Output:
% data_epoch_spectra is channels_cmb X epochs X frequencies for 'sparse'(default)
%   or
% data_epoch_spectra is channels X channels X epochs X frequencies for 'full'
%
% Example:
% % Output sparse matrix without auto-spectra
% [f,corss_data_epoch_spectra,data_epoch] = ...
%     ecog_crossspectra(data_epoch,stims,fft_w,fft_t,fft_ov,srate,regress_erp,f_limits,'sparse')
% % Output full matrix with auto-spectra
% [f,corss_data_epoch_spectra,data_epoch] = ...
%     ecog_crossspectra(data_epoch,stims,fft_w,fft_t,fft_ov,srate,regress_erp,f_limits,'full')
% % Output sparse matrix with auto-spectra
% [f,corss_data_epoch_spectra,data_epoch] = ...
%     ecog_crossspectra(data_epoch,stims,fft_w,fft_t,fft_ov,srate,regress_erp,f_limits,'sparse+')
% % Output magnitude-squared coherence isntead of cross-spectrum
% [f,corss_data_epoch_spectra,data_epoch] = ...
%     ecog_crossspectra(data_epoch,stims,...,'coh')
%
%
% K.Yuasa, 20201116 - modify from ecog_spectra
% K.Yuasa, 20210511 - enable to set limit in output frequency
% K.Yuasa, 20210514 - change output data structure
% K.Yuasa, 20230315 - add coherence mode

% regress erp out
if ~exist('f_lim','var')||isempty(f_lim)
    f_lim = [0 inf];
end
if ~exist('outtype','var')||isempty(outtype)
    outtype = 'sparse';
end
if ~exist('calcmode','var')||isempty(calcmode)
    calcmode = 'default';
end
if reg_erp==1 
    data_epoch_orig = data_epoch;
    for k = 1:size(data_epoch,1)%channels
        disp(['regress erp el ' int2str(k)])
        for m = 1:size(data_epoch,2)%epochs
            x = squeeze(data_epoch(k,m,:));
            if all(isnan(x))
                % regress ERP out
                s = stims(m);
                av_erp = squeeze(mean(data_epoch_orig(k,stims==s,:),2,'omitnan'));
                [~,~,reg_R] = regress(x,av_erp);
            else
                reg_R = x;
            end
            data_epoch(k,m,:) = reg_R;
        end
    end
end

% calculate spectra to get length of f to initialize spectra
[ix,iy]=find(~any(isnan(data_epoch(:,:,fft_t)),3),1);       % 1st non-nan indices
[~,f] = pwelch(squeeze(data_epoch(ix,iy,fft_t)),fft_w,fft_ov,srate,srate);
fidx = f >= f_lim(1) & f <= f_lim(2);
f = f(fidx);
switch outtype
    case 'sparse',        shifttri = 1;
    case 'sparse+',       shifttri = 0;
end
switch outtype
    case {'sparse','sparse+'}
        [chanx,chany]=meshgrid(1:size(data_epoch,1),1:size(data_epoch,1));
        chanl = true(size(chanx));
        chancmb = [chanx(tril(chanl,-shifttri)), chany(tril(chanl,-shifttri))];
        data_epoch_spectra = zeros(size(chancmb,1),size(data_epoch,2),length(f));
    case 'full'
        data_epoch_spectra = zeros(size(data_epoch,1),size(data_epoch,1),size(data_epoch,2),length(f));
    otherwise
        error('%s is unknown input',outtype);
end
switch lower(calcmode)
    case {'default','cpsd'},         func_cross = @cpsd;     func_auto = @pwelch; 
    case {'coh','mscoh','mscohere'}, func_cross = @mscohere; func_auto = @(varargin) mscohere(varargin{1},varargin{1},varargin{2:end}); 
end

% calculate cross-spectra (skip k==l)
ii=1;
for k = 1:size(data_epoch,1)%channels
    disp(['fft el ' int2str(k) ' x other']);
    switch outtype
    case {'sparse','sparse+'}
            for l = (k+shifttri):size(data_epoch,1)%channels
                x = permute(data_epoch(k,:,:),[3,2,1]);
                y = permute(data_epoch(l,:,:),[3,2,1]);
                nanepoch = or(all(isnan(x),1),all(isnan(y),1));
                [Pxy] = func_cross(x(fft_t,~nanepoch),y(fft_t,~nanepoch),fft_w,fft_ov,srate,srate);
                data_epoch_spectra(ii,~nanepoch,:) = Pxy(fidx,:)';
                ii = ii + 1;
            end
        case 'full'
            for l = 1:size(data_epoch,1)%channels
                if l<k
                    data_epoch_spectra(k,l,:,:) = conj(data_epoch_spectra(l,k,:,:));
                elseif l>k
                    x = permute(data_epoch(k,:,:),[3,2,1]);
                    y = permute(data_epoch(l,:,:),[3,2,1]);
                    nanepoch = or(all(isnan(x),1),all(isnan(y),1));
                    [Pxy] = func_cross(x(fft_t,~nanepoch),y(fft_t,~nanepoch),fft_w,fft_ov,srate,srate);
                    data_epoch_spectra(k,l,~nanepoch,:) = Pxy(fidx,:)';
                elseif l==k
                    x = permute(data_epoch(k,:,:),[3,2,1]);
                    nanepoch = all(isnan(x),1);
                    [Pxx] = func_auto(x(fft_t,~nanepoch),fft_w,fft_ov,srate,srate);
                    data_epoch_spectra(k,k,~nanepoch,:) = Pxx(fidx,:)';
                end
            end
    end
end
