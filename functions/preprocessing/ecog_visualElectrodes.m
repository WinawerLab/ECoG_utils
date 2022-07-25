function [channels,chanidx] = ecog_visualElectrodes(channels,visualfields,checknode)
% [channels,index] = ecog_visualElectrodes(channels,[visualfields, checknode])
%   outputs channels which are associated to any visual areas
% 
% Input
%     channels:  	 channel table
%     visualfields:  field names to recognize visual areas
%     checknode:     true, false(defalut)
%                    if true, check 'matchednode' field in channels at first
% 
% Output
%     channels:  	channel table
%     index:        logical array to match input and output channels
% 
% See also, ecog_HDgridElectrodes
% 
% K.Yuasa, 2022

% Atlas fileds
if ~exist('checknode','var') || isempty(checknode)
    checknode = false;
end
if ~exist('visualfields','var') || isempty(visualfields)
    checknode = true;
    visualfields = {'bensonarea','wangarea','hcparea',...
                    '*benson14_varea','*wang15_mplbl','benson20_mplbl',...
                    'wangprob_*','hcpprob_*',...
                    '*wang15_fplbl_*','benson20_fplbl_*'};
end
visualfields = regexptranslate('wildcard',string(visualfields));

% Check node
chanidx    = [];
chfldnames = reshape(fieldnames(summary(channels)),1,[]);
if checknode
    
    curfld = "matchednode";
    fldlist = ismember(curfld,chfldnames);
    if fldlist
        chanidx = ~isnan(channels.(curfld));
    end
    
end

% Check atlas
if isempty(chanidx)
    
isnumflds  = all(cellfun(@isnumeric,table2cell(channels)),1);
for curfld = reshape(visualfields,1,[])
    fldlist = ~cellfun(@isempty,regexp(chfldnames,curfld));
    if any(fldlist&isnumflds)
        visch  = sum(channels{:,fldlist&isnumflds},2)>0;
    elseif any(fldlist&~isnumflds)
        visch  = any(~contains(channels{:,fldlist&~isnumflds},'none'),2);
    else
        visch  = false(height(channels),1);
    end
    chanidx = [chanidx, visch];
end

end

% Output
chanidx  = any(chanidx,2);
channels = channels(chanidx,:);

