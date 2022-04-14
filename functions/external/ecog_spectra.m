function [f,data_epoch_spectra,data_epoch] = ...
    ecog_spectra(data_epoch,stims,fft_w,fft_t,fft_ov,srate,reg_erp)

% Function to calculate powerspectrum for across several ECoG channels and
% epochs using Matlab's pwelch function.
% 
% Input:
% data_epoch: channels X epochs X time
%
% stims: if regress_erp==1, the ERP calculated per condition can be regressed out.
%   The condition is indicated in stims, a vector with one value for each
%   epoch. If choosing to regress erp out, make sure data are baseline
%   corrected first, so the ERP can be calculated correctly!
%
% Output:
% data_epoch_spectra is channels X epochs X frequencies
%
% Example:
% [f,data_epoch_spectra,data_epoch] = ...
%     ecog_spectra(data_epoch,stims,fft_w,fft_t,fft_ov,srate,regress_erp)
%
%
% Dora Hermes, 2017

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
[~,f] = pwelch(squeeze(data_epoch(1,1,fft_t)),fft_w,fft_ov,srate,srate);
data_epoch_spectra = zeros(size(data_epoch,1),size(data_epoch,2),length(f));
clear Pxx

% calculate powerspectra
for k = 1:size(data_epoch,1)%channels
    disp(['fft el ' int2str(k)])
    x = permute(data_epoch(k,:,:),[3,2,1]);     
    nanepoch = all(isnan(x),1);
    [Pxx,f] = pwelch(x(fft_t,~nanepoch),fft_w,fft_ov,srate,srate);
    data_epoch_spectra(k,~nanepoch,:) = Pxx';
end