function data = ecog_notch(data,srate,notch_freq,nharmonic,norder)

% data = ecog_notch(data,srate,notch_freq,[harmonics,order])
%   apply notch filter around notch_freq Hz and its harmonics.
% 
% INPUTS:
%   data:       ( time X channels )    
%   srate:      sampling rate (Hz)
%   notch_freq: notch frequency (Hz) [default = 60]
%   harmonic:   Nth-harmonics to be notched [default = 3]
%   order:      Nth-order filter [default = 5]

% 20220824 Yuasa

%-- Set parameters
narginchk(2,inf);
if ~exist('notch_freq', 'var') || isempty(notch_freq), notch_freq = 60; end   % notch frequency
if ~exist('nharmonic', 'var') || isempty(nharmonic), nharmonic = 3; end   % Nth-harmonic
if ~exist('norder', 'var') || isempty(norder), norder = 5; end   % Nth-order

nharmonic = round(nharmonic(1));
norder    = round(norder(1));

%-- Expand notch frequencies to include all harmonics
notch_freq = reshape(notch_freq(:).*(1:nharmonic),1,[]);

%-- Apply notch filter
fprintf('[%s] Notching out ', mfilename);
for inotch = notch_freq
    fprintf('%g ', inotch);
    [ib, ia]=butter(norder,([-1 1]+inotch)./srate*2,'stop'); %60hz
    data=filtfilt(ib,ia,data);
end
fprintf('Hz\n');
