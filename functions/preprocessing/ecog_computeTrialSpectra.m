function [f,trials_spectra] = ecog_computeTrialSpectra(trials, specs)

% Calculate powerspectrum for all trials and channels in 'trials' struct
% using Dora's script which uses Matlab's pwelch function.
% 
% Input:
% trials: structure as outputted by preprocess_bairECOG
% specs: specification of fft parameters for ecog_spectra.m function

% FFT SETTINGS
fft_w = window(@hann,specs.window); % window width
fft_ov = specs.ov; % overlap
reg_erp = specs.reg_erp;

% TRIALS
data = trials.evoked;
data_epoch = zeros(size(data,1), size(data,3), size(data,2));
for ii = 1:size(data,3)
    data_epoch(:,ii,:) = data(:,:,ii);
end

% TIME
stims = ones(1,size(data_epoch,2));
fft_t = trials.time > specs.t(1) & trials.time < specs.t(2); 
srate = trials.fsample;

% COMPUTE
[f,trials_spectra] = ecog_spectra(data_epoch,stims,fft_w,fft_t,fft_ov,srate,reg_erp);

end