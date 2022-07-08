function [bb_amp_low,alpha_amp,alpha_freq,alpha_width,fit_fd2,out_exp,beta_params] = ...
    ecog_fitalpha(f,f_use4fit,f_alpha,data_base,data_fit,neg_amp,fit_beta)

% function fits broadband + gaussian for <data_base - data_fit>
% [bb_amp_low,alpha_amp,alpha_freq,alpha_width,fit_fd2,out_exp,[beta_params]] = ...
%     ecog_fitalpha(f,f_use4fit,[f_alpha],data_base,data_fit,[neg_amp,est_beta]);
%
% input:
% f: frequencies
% f_use4fit: frequencies to use fitting (ex. [1:30])
% f_alpha: boundary for peak alpha frequnecy (default: 8-13Hz allowing -2,+4Hz extension)
% data_base: used to fit exp: (1/f^exp) - enter power (not log)
% data_fit: used to fit weights and gaussian - enter power (not log)
% neg_amp: if false (default), only fit to alpha suppression
%          if true, fit both to alpha suppression and enhancement
% fit_beta:if false (default), never consider beta bump
%          if true, fit beta bump in 14-30Hz
%
% output (weight_pwr_at_alpha_freq weight_gauss alpha_freq width_gauss fit_fd2 exp)

% 20190903 Yuasa: modified from ecog_fitgamma.m in gammaModel toolbox (Hermes, 2015)
% 20191112 Yuasa: update parameters & do fitting twice
% 20191114 Yuasa: loop second fitting & introduce neg_amp
% 20191115 Yuasa: introduce multiple initial points for alpha fitting
%                 introduce f_alpha
% 20200224 Yuasa: simplify the script
% 20210907 Yuasa: enable to fit beta

if nargin<7,    fit_beta = false;    end
if nargin<6,    neg_amp = false;    end
if nargin<5
    data_fit    = data_base;
    data_base   = f_alpha;
    f_alpha     = [];
elseif numel(data_base) ~= numel(data_fit)
    fit_beta    = neg_amp;
    neg_amp     = data_fit;
    data_fit    = data_base;
    data_base   = f_alpha;
    f_alpha     = [];
end
f           = f(:);
data_fit    = data_fit(:);
data_base   = data_base(:);
f_sel = ismember(f,f_use4fit);
x_in = (log10(data_base(f_sel)) - log10(data_fit(f_sel)));
f_in = log10(f(f_sel));
func_model = @(X,F,C,X2) X(1) - C*F + X(2)*sqrt(2*pi)*normpdf(F,X(3),X(4)) + X2(1)*sqrt(2*pi)*normpdf(F,X2(2),X2(3));
%%% func_powerBump(X,x_in,f_in,x_slope) == x_in - func_model(X,f_in,x_slope)

if isempty(f_alpha)
    f_bounds    = [6 17];
    f_cores     = [8 13];
elseif numel(f_alpha)==1
    f_bounds    = [f_alpha f_alpha];
    f_cores     = f_bounds;
elseif isvector(f_alpha)
    f_bounds    = f_alpha;
    f_cores     = f_bounds;
else
    f_bounds    = f_alpha(1,:);
    f_cores     = f_alpha(2,:);
    assert(f_cores(1)>=f_bounds(1) && f_cores(2)<=f_bounds(2), 'f_alpha is invalid');
end
if neg_amp,     limminamp = -Inf;
else,           limminamp = 0;
end

%-- fit exponent to base
x_slope = 0;
f_base = f_sel & ~((f>=6 & f<=15) | (f>20 & f<30));
p = polyfit(log10(f(f_base)),(log10(data_base(f_base)) - log10(data_fit(f_base))),1);
x_slope = [x_slope -p(1)];
f_base = f_sel & ~(f>=6);
p = polyfit(log10(f(f_base)),(log10(data_base(f_base)) - log10(data_fit(f_base))),1);
x_slope = [x_slope -p(1)];

%-- cost function (set priority for f_cores) 2≥cost≥1
maxdist = max(abs(log10(f_bounds) - log10(f_cores)));
if ~maxdist, maxdist = 1; end
func_cost = @(F) double(prod(F - log10(f_cores))>0).*min(abs(F - log10(f_cores)))./maxdist+1;

%-- fit powerlaw and gaussian and plot
my_options = optimset('Display','off','Algorithm','trust-region-reflective');
if fit_beta
nparams = 8;
x0 = cell(1,nparams);
[x0{:}] = ndgrid([0],[0],log10([8 10 13]),log10([15/8 15/10 15/13])./4,x_slope,[0],log10([20 30]),log10([30/20])./4);
x0      = reshape(cat(nparams+1,x0{:}),[],nparams);
[x,fiterr] = mlsqnonlin(@(X) func_cost(X(3)).*(x_in - func_model(X(1:4),f_in,X(5),X(6:8))),x0,...
    [-Inf limminamp log10(f_bounds(1)) (log10(13/12))/4 -Inf 0 log10(max(15,f_bounds(2)+2)) (log10(31/30))/4],...
    [Inf Inf log10(f_bounds(2)) (log10(17/6))/4 Inf Inf log10(max(30,f_bounds(2)+4)) (log10(34/15))/4],...
    my_options);
else
nparams = 5;
x0 = cell(1,nparams);
[x0{:}] = ndgrid([0],[0],log10([8 10 13]),log10([15/8 15/10 15/13])./4,x_slope);
x0      = reshape(cat(nparams+1,x0{:}),[],nparams);
[x,fiterr] = mlsqnonlin(@(X) func_cost(X(3)).*(x_in - func_model(X(1:4),f_in,X(5),[0,0,1])),x0,...
    [-Inf limminamp log10(f_bounds(1)) (log10(13/12))/4 -Inf],...
    [Inf Inf log10(f_bounds(2)) (log10(17/6))/4 Inf],...
    my_options);
end

x_slope     = x(5);
alpha_amp   = -x(2);
alpha_freq  = x(3);
alpha_width = x(4);
bb_amp_low  = -(x(1)-x_slope*alpha_freq);

out_exp = [x_slope x(1)];

%-- fit to pre-post data in log-space
% fit_fd2 = x(1) - x_slope*log10(f) + ...
%     -alpha_amp*sqrt(2*pi) * normpdf(log10(f),alpha_freq,alpha_width);
if fit_beta
    beta_params = x(6:8);
fit_fd2 = [func_model(x(1:4),log10(f),x_slope,beta_params),...
           func_model(x(1:4),log10(f),x_slope,[0,0,1])];
else
    beta_params = nan(1,3);
fit_fd2 = func_model(x(1:4),log10(f),x_slope,[0,0,1]);
end

function  [xCurrent,Resnorm,FVAL,EXITFLAG,OUTPUT,LAMBDA,JACOB] = mlsqnonlin(FUN,xCurrents,LB,UB,options,varargin)
% MLSQNONLIN computes LSQNONLIN with multiple initial points, and returns
% the outputs with the best Resnorm.
% 
% xCurrents: MxN matrix with M initial points with N parameters
% 
% See also: LSQNONLIN

outs = cell(1,7);
minResnorm = Inf;
for initp=1:size(xCurrents,1)
    [outs{:}] = lsqnonlin(FUN,xCurrents(initp,:),LB,UB,options,varargin{:});
%     disp(outs{2});
    if outs{2} < minResnorm
        [xCurrent,Resnorm,FVAL,EXITFLAG,OUTPUT,LAMBDA,JACOB] = outs{:};
        minResnorm = Resnorm;
%         disp([xCurrents(initp,:);xCurrent]);
    end
end

