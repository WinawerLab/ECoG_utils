function [out_stat,upci_base] = ecog_bootstrapStat(data_bs,base_bs,option_test)
% Function to calculate the p-value (out_stat) for the difference between
% two bootstrapped distributions.
%
% Input 
% data_bs: stimulus bootstraps (stimuli X number of bootstraps)
% base_bs: baseline bootstraps (1 X number of bootstraps)
% option_test: 'bootdiff' for getting p values based on the difference
% between every bootstrap and the baseline
%
% Output
% out_stat depends on option_test:
%   bootdiff: p values
%   zstat: t statistic
% 
% Dora Hermes, 2019

out_stat = zeros(size(data_bs,1),1);

if isequal(option_test,'zstat')
    
    for kk = 1:size(data_bs,1) % number of stimuli
        % these are means of each bootstrap: 
        xmn = data_bs(kk,:);
        ymn = base_bs;

        % now we want to compute statistics
        mn1 = mean(xmn);
        mn2 = mean(ymn);
        se1 = std(xmn);  % this is already like a "standard error"
        se2 = std(ymn);
        out_stat(kk) = (mn1-mn2)./sqrt(se1.^2+se2.^2); % take the standard errors and combine them
    end
    
elseif isequal(option_test,'bootdiff')
    
	base_ci = prctile(base_bs,[16 84]);
    for kk = 1:size(data_bs,1) % number of stimuli
        boot_stat = data_bs(kk,:)-base_bs;
        out_stat(kk) = sum(boot_stat<0)/size(data_bs,2);
    end
    upci_base = base_ci(2);
end
