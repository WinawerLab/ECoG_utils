function [out_exp,bb_amp,gamma_amp,gamma_freq,gamma_width,fit_f2] = ...
    ecog_fitgamma(f,f_use4fit,data_base,data_fit)

% function fits broadband + gaussian
% [out_exp,bb_amp,gamma_amp,gamma_freq,gamma_width,fit_f2] = ...
%     ecog_fitgamma(f,f_use4fit,data_base,data_fit);
%
% input:
% f: frequencies
% f_use4fitb= [35:57 65:115 126:175 186:200];
% data_base: used to fit exp: (1/f^exp) - enter power (not log)
% data_fit: used to fit weights and gaussian - enter power (not log)
%
% output (exp weight_pwr weight_gauss gamma_freq fit_f2)
% D Hermes, 2015

f_sel = ismember(f,f_use4fit);
x_base = data_base(f_sel);
x_in = data_fit(f_sel);
f_in = f(f_sel);

% fit exponent to base
p = polyfit(log10(f_in),log10(x_base)',1);
out_exp = [-p(1) p(2)];

% fit powerlaw and gaussian and plot
my_options = optimset('Display','off','Algorithm','trust-region-reflective');
[x] = lsqnonlin(@(x) func_powerBump(x,log10(x_in),log10(f_in'),out_exp(1)),...
    [0 0 log10(40) .05],[-Inf 0 log10(30) .03],[Inf Inf log10(80) .08],...
    my_options);
bb_amp = x(1);
gamma_amp = x(2);
gamma_freq = x(3);
gamma_width = x(4);

% fit to data in log-space
fit_f2 = bb_amp-out_exp(1)*log10(f) + ...
    gamma_amp*sqrt(2*pi) * normpdf(log10(f),gamma_freq,gamma_width);

%%
% 
% figure,hold on
% plot(f_in,log10(x_base),'k.','LineWidth',2)
% plot(f_in,log10(x_in),'r.','LineWidth',2)
% plot(f_in,fit_f2(f_sel),'-','Color',[.5 .5 .5])
% 
% figure,hold on
% plot(log10(f_in),log10(x_base),'k.','LineWidth',2)
% plot(log10(f_in),log10(x_in),'r.','LineWidth',2)
% plot(log10(f_in),fit_f2(f_sel),'-','Color',[.5 .5 .5])