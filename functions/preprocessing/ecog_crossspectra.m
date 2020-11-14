function [f,data_epoch_spectra,data_epoch] = ...
    ecog_crossspectra(data_epoch,stims,fft_w,fft_t,fft_ov,srate,reg_erp)

% Function to calculate cross-powerspectrum for across several ECoG
% channels and epochs using Matlab's cpsd function.
% 
% Input:
% data_epoch: channels X epochs X time
%
% stims: if repress_erp==1, the ERP calculated per condition can be regressed out.
%   The condition is indicated in stims, a vector with one value for each
%   epoch. If choosing to regress erp out, make sure data are baseline
%   corrected first, so the ERP can be calculated correctly!
%
% Output:
% data_epoch_spectra is channels X channels X epochs X frequencies
%
% Example:
% [f,corss_data_epoch_spectra,data_epoch] = ...
%     ecog_crossspectra(data_epoch,stims,fft_w,fft_t,fft_ov,srate,reg_erp)
%
%
% Ken Yuasa, 2020 (based on ecog_spectra)

% regress erp out
if reg_erp==1 
    data_epoch_orig = data_epoch;
    for k = 1:size(data_epoch,1)%channels
        disp(['regress erp el ' int2str(k)])
        for m = 1:size(data_epoch,2)%epochs
            x = squeeze(data_epoch(k,m,:));
            if all(isnan(x))
                % regress ERP out
                s = stims(m);
                av_erp = squeeze(nanmean(data_epoch_orig(k,stims==s,:),2));
                [~,~,reg_R] = regress(x,av_erp);
            else
                reg_R = x;
            end
            data_epoch(k,m,:) = reg_R;
        end
    end
end

% calculate spectra to get length of f to initialize spectra
[~,f] = pwelch(squeeze(data_epoch(1,1,fft_t)),fft_w,fft_ov,srate,srate);
data_epoch_spectra = zeros(size(data_epoch,1),size(data_epoch,1),size(data_epoch,2),length(f));

% calculate cross-spectra (skip k==l)
for k = 1:size(data_epoch,1)%channels
    disp(['fft el ' int2str(k) ' x other']);
    for l = 1:size(data_epoch,1)%channels
        if l<k
            data_epoch_spectra(k,l,:,:) = conj(data_epoch_spectra(l,k,:,:));
        elseif l>k
            x = squeeze(data_epoch(k,:,:))';
            y = squeeze(data_epoch(l,:,:))';
            nanepoch = or(all(isnan(x),1),all(isnan(y),1));
            [Pxy,f] = cpsd(x(fft_t,~nanepoch),y(fft_t,~nanepoch),fft_w,fft_ov,srate,srate);
            data_epoch_spectra(k,l,~nanepoch,:) = Pxy';
        end
    end
end