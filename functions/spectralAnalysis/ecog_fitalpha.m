function [bb_amp_low,alpha_amp,alpha_freq,alpha_width,fit_fd2,out_exp,beta_params,fiterr,gof] = ...
    ecog_fitalpha(f,f_use4fit,f_alpha,data_base,data_fit,neg_amp,fit_beta)

% function fits broadband + gaussian for <data_base - data_fit>
% [bb_amp_low,alpha_amp,alpha_freq,alpha_width,fit_fd2,out_exp,[beta_params,fiterr,cod]] = ...
%     ecog_fitalpha(f,f_use4fit,[f_alpha],data_base,data_fit,[neg_amp],[est_beta]);
%
% input:
% f: frequencies
% f_use4fit: frequencies to use fitting (ex. [1:30])
% f_alpha: boundary for peak alpha frequnecy (default: 8-13Hz allowing -2,+4.5Hz extension)
% data_base: used to fit exp: (1/f^exp) - enter power (not log)
% data_fit: used to fit weights and gaussian - enter power (not log)
% neg_amp: if false (default), only fit to alpha suppression
%          if true, fit both to alpha suppression and enhancement
% fit_beta:if false (default), never consider beta bump
%          if true, fit beta bump in 15-30Hz
%          if -1,   fit theta bump in 3–6Hz instead of beta
%
% output (weight_pwr_at_alpha_freq weight_gauss alpha_freq width_gauss fit_fd2 exp)

% 20190903 Yuasa: modified from ecog_fitgamma.m in gammaModel toolbox (Hermes, 2015)
% 20191112 Yuasa: update parameters & do fitting twice
% 20191114 Yuasa: loop second fitting & introduce neg_amp
% 20191115 Yuasa: introduce multiple initial points for alpha fitting
%                 introduce f_alpha
% 20200224 Yuasa: simplify the script
% 20210907 Yuasa: enable to fit beta
% 20220829 Yuasa: enable to compute goodness of fitting
% 20220830 Yuasa: modify parameters
% 20220905 Yuasa: add theta fit option

if nargin<7,    fit_beta = false;    end
if nargin<6,    neg_amp = false;    end
if nargin<5
    data_fit    = data_base;
    data_base   = f_alpha;
    f_alpha     = [];
elseif numel(data_base) ~= numel(data_fit) && numel(data_base) == numel(f_alpha)
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

if isempty(f_alpha)     % default = [8 13]
    f_bounds    = [6 17.5];
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

%-- fit exponent to base
x_slope = 0;    x_intcp = min(x_in);
f_base = f_sel & ~((f>=6 & f<=15) | (f>20 & f<30));
p = polyfit(log10(f(f_base)),(log10(data_base(f_base)) - log10(data_fit(f_base))),1);
x_slope = [x_slope -p(1)];  x_intcp = [x_intcp p(2)];
f_base = f_sel & ~(f>=6);
p = polyfit(log10(f(f_base)),(log10(data_base(f_base)) - log10(data_fit(f_base))),1);
x_slope = [x_slope -p(1)];  x_intcp = [x_intcp p(2)];

%-- cost function (set priority for f_cores) 2≥cost≥1
maxdist = max(abs(log10(f_bounds) - log10(f_cores)));
if ~maxdist, maxdist = 1; end
func_cost  = @(F) 1 + double(prod(F - log10(f_cores))>0).*min(abs(F - log10(f_cores)))./maxdist;
%-- cost function for beta fitting: fit with a bump is better
func_costb = @(F1,F2,A1,A2) 1 + double(A2>abs(A1./10) & F2<(F1+log10(19/13)) & F2>(F1-log10(10/5)))...
                              + double(A2>abs(A1./5) & F2<(F1+log10(16/13)) & F2>(F1-log10(8/6)))...
                              + double(A2>abs(A1.*0.3));
%-- restriction for beta fitting: apply zero if F2 is not greater than F1 or A1 is at least a half of A2
func_restr = @(F1,F2,A1,A2) double((F2>F1 | (F2<F1 & fit_beta<0)) & abs(2*A1)>A2);
%-- fitting function
func_fit   = @(X,F,data) func_costb(X(3),X(7),X(2)./X(4),X(6)./X(8)) .* func_cost(X(3)).*...
    (data - func_model(X(1:4),F,X(5),X(6:8)).*func_restr(X(3),X(7),X(2)./X(4),X(6)./X(8)));

%-- fit powerlaw and gaussian and plot
my_options = optimset('Display','off','Algorithm','trust-region-reflective');
nparams = 5;
x0     = setinitval(nparams,x_slope,x_intcp);
xlimit = setfitlim(nparams,f_bounds,neg_amp);
[x,fiterr] = mlsqnonlin(@(X) func_fit([X [0,Inf,1]],f_in,x_in),x0,...
    xlimit(1,:),xlimit(2,:),...
    my_options);
gof = calcgof(fiterr,x_in,f_in);    % this gof include cost function

if fit_beta
nparams = 8;
x0     = [x [0,Inf,1]];
if gof<20, x = []; end
x0     = [x0;setinitval(nparams,x_slope,x_intcp,x,fit_beta==-1)];    % change mode to theta if fit_beta == -1
xlimit = setfitlim(nparams,f_bounds,neg_amp,fit_beta==-1);
[x,fiterr] = mlsqnonlin(@(X) func_fit(X,f_in,x_in),x0,...
    xlimit(1,:),xlimit(2,:),...
    my_options);
end

alpha_amp   = -x(2);
alpha_freq  = x(3);
alpha_width = x(4);
bb_amp_low  = -(x(1)-x(5)*alpha_freq);

out_exp = x([5,1]);

%-- fit to pre-post data in log-space
if nargout > 6
if fit_beta
    beta_params = x(6:8);
    fit_fd2 = [func_model(x(1:4),log10(f),x(5),beta_params),...
               func_model(x(1:4),log10(f),x(5),[0,0,1])];
else
    beta_params = nan(1,3);
    fit_fd2 = func_model(x(1:4),log10(f),x(5),[0,0,1]);
end
end

%-- compute goodness of fitting
if nargout > 8
    gof = calcgof(fit_fd2(f_sel,1),x_in,f_in);
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

function x0 = setinitval(nparams,x_slope,x_intcp,x0pre,fit_theta)
% Set initial values for fitting
if ~exist('fit_theta','var') || isempty(fit_theta)
    fit_theta = false;
end

nbase = length(x_slope);
assert(length(x_intcp)==nbase,'Sizeds of x_slope and x_intcp must be the same');
if ~exist('x0pre','var'), x0pre =[];  end
nparamspre = size(x0pre,2);
nnewparams = nparams-nparamspre;

x0 = cell(1,nnewparams);
if fit_theta        % fit theta instead of beta
parameterlist = {[0],[0.1],log10([8 10 13]),log10([15/8 15/10 15/13])./4,[1:nbase],[0 0.05],log10([3 5]),log10([8/3])./4};
else                % fit beta
parameterlist = {[0],[0.1],log10([8 10 13]),log10([15/8 15/10 15/13])./4,[1:nbase],[0 0.05],log10([20 30]),log10([30/20])./4};
end
[x0{:}] = ndgrid(parameterlist{(nparamspre+1):nparams});
x0      = reshape(cat(nnewparams+1,x0{:}),[],nnewparams);
x0      = horzcat(repmat(x0pre,size(x0,1),1), x0);
if isempty(x0pre)
x0(:,2) = x0(:,2).*x0(:,4);  % correct alpha amplitude
x0(:,[5,1])  = [x_slope(x0(:,5))', x_intcp(x0(:,5))'];  % correct baseline properties
end
if nparams>5
x0(:,6) = x0(:,6).*x0(:,8);  % correct beta amplitude
end

function xlimit = setfitlim(nparams,f_bounds,neg_amp,fit_theta)
% Set limit values for fitting
if ~exist('fit_theta','var') || isempty(fit_theta)
    fit_theta = false;
end
if ~exist('neg_amp','var') || isempty(neg_amp)
    neg_amp   = false;
end

if neg_amp,     limminamp = -Inf;
else,           limminamp = 0;
end

if fit_theta        % fit theta instead of beta
    xlimit = [...
            [-Inf limminamp log10(f_bounds(1)) (log10(14/12))/4 -Inf 0   log10(min(3,f_bounds(1)-2)) (log10(8/6))/4];
            [Inf  Inf       log10(f_bounds(2)) (log10(18/4))/4   Inf Inf log10(6)                    (log10(10/2))/4]...
            ];
else                % fit beta
    xlimit = [...
            [-Inf limminamp log10(f_bounds(1)) (log10(14/12))/4 -Inf 0   log10(15)                    (log10(32/30))/4];
            [Inf  Inf       log10(f_bounds(2)) (log10(18/4))/4   Inf Inf log10(max(30,f_bounds(2)+2)) (log10(35/10))/4]...
            ];
end
xlimit = xlimit(:,1:nparams);


function gof = calcgof(fiterr,data,f)
% Compute goodness of fitting
%   gof = calcgof(fiterr,data,f)
%   gof = calcgof(model,data,f)

if isequal(size(fiterr),size(data))
    % interpret 1st input as model instead of fiterr
    fiterr = sum((data-fiterr).^2,'omitnan');
end
    
p = polyfit(f,data,1);
base = p(1).*f+p(2);
gof = 100*(1 - fiterr ./sum((data-base).^2,'omitnan'));
    