function [coefs,predictor] = ecog_regressERP(data_epoch,t_base,stims,stim4predict,maxlag,negcoef4lag,t_ref)
% function to regress the ERP from some data
% USAGE:
% [coef,predictor] = ecog_regressERP(data_epoch,t_base,stims,stim4predict)
% [coef,predictor] = ecog_regressERP(data_epoch,t_base,stims,stim4predict,maxlag,negcoef4lag,[t_ref])
% 
% Inputs: 
%   data_epoch   % electrodes X time X epochs
%   t_base       % logical array to indicate the baseline time points
%   stims        % different code for different conditions
%   stim4predict % the number of the condition code to construct the predictor
%   maxlag       % allows lag in the range from -maxlag to maxlag for regression
%                  if 0, regression is applied without lag (default: [0])
%                  (ex: [5] for maximum 10 ms lag at 500 Hz sampling rate)
%   negcoef4lag  % if false, adopt the lag when the coefficient is largest (default)
%                  if true, adopt the lag when the absolute coefficient is largest
%   t_ref        % logical array to indicate the reference time points to apply regression (default: all time points)
% 
% Outputs: 
%   coef         % coefficients of regressor
%   predictor    % regressor (ERP)
% 
% See also, ecog_regressout

% 20190903 Yuasa: modified from ecog_spectra.m in ECoG_utils
% 20200303 Yuasa: make computations effective & use 'parfor'
% 20210413 Yuasa: use same predictor for all stims
% 20210425 Yuasa: add maxlag option
% 20220825 Yuasa: alias to ecog_regressout

%-- check data_epoch
narginchk(4,6);
if nargin < 7,  t_ref = [];             end
if nargin < 6,  negcoef4lag = false;    end
if nargin < 5,  maxlag = 0;             end

[~,coefs,predictor] = ecog_regressout(data_epoch,t_base,stims,maxlag,t_ref,negcoef4lag,stim4predict);