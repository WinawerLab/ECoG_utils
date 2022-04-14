function [atlasFile, isnorm, normthresh, issmry, atlasName] = interpretAtlasNames(atlasName,flg)

% [atlasFile, isnorm, normthresh, issmry, atlasName] = INTERPRETATLASNAMES(atlasName)
%   Interpret atlas names with normalize options.
% Example 1: 
%   INTERPRETATLASNAMES(smry_wang15_fplbl_norm_5)
%   atlasName = 'smry_wang15_fplbl_norm_5' is interepreted as 
%       atlasFile   = 'wang15_fplbl'
%       isnorm      = true;
%       normthresh  = 5;
%       issmry      = true;
%       atlasName   = 'smry_wang15_fplbl_norm_5'
% 
% Example 2: 
%   INTERPRETATLASNAMES(wang15_fplbl_norm)
%   atlasName = 'wang15_fplbl_norm' is interepreted as 
%       atlasFile   = 'wang15_fplbl'
%       isnorm      = true;
%       normthresh  = 0;
%       issmry      = false;
%       atlasName   = 'wang15_fplbl_norm_0'

% K.Yuasa, BAIR 2022

cellflg = iscell(atlasName);
if ~cellflg, atlasName = {atlasName};  end
atlasNum       = numel(atlasName);
atlasFile      = cell(1,atlasNum);
isnorm         = false(1,atlasNum);
normthresh     = zeros(1,atlasNum);
normthreshSet  = cell(1,atlasNum);
issmry         = false(1,atlasNum);

if nargin > 1 && ischar(flg) && strcmpi(flg,'rename')
    dorename = true;
else
    dorename = false;
end

for a = 1:atlasNum
    % Look up area names
    currentAtlas = atlasName{a};
    
    % Check normalized
    normchk = contains(currentAtlas,'_fplbl_norm');
    if normchk
        [atlastok] = regexp(currentAtlas,'_norm_?','split');
        baseName   = atlastok{1};
        atlastok   = strsplit(strjoin(atlastok(2:end),'_'),'_');
        % Check normalized threshold
        thresh   = cellfun(@str2double,atlastok);
        thresh(isnan(thresh)) = 0;
        thresh = reshape(unique(thresh),1,[]);
        currentAtlas = regexp(currentAtlas,'.*_fplbl_norm','match','once');
    else
        baseName   = currentAtlas;
        thresh  = 0;
    end
    % Check summary label
    smrychk  = contains(currentAtlas,'smry_');
    baseName = strrep(baseName, 'smry_', '');
    
    % Rename
    if dorename
        currentAtlas = regexprep(currentAtlas,'_fplbl_norm.*','_mplbl');
        currentAtlas = strrep(currentAtlas,'smry_','');
    end
    
    % Outputs
    atlasName{a}       = currentAtlas;
    atlasFile{a}       = baseName;
    isnorm(a)          = normchk;
    normthreshSet{a}   = thresh;
    issmry(a)          = smrychk;
end

if cellflg
% Expand normthresh unless input is not cell
    for a = atlasNum:-1:1
        Nthresh = length(normthreshSet{a});
        % Check necessity of expansion
        if Nthresh > 1
            if ~dorename && isnorm(a)
            TMPatlasFile   = repmat(atlasFile(a),1,Nthresh);
            TMPisnorm      = repmat(isnorm(a),1,Nthresh);
            TMPnormthresh  = normthreshSet{a};
            TMPissmry      = repmat(issmry(a),1,Nthresh);
            TMPatlasName   = arrayfun(@(d) sprintf('%s_%d',atlasName{a},round(d)),...
                                      normthreshSet{a},'UniformOutput',false);
            else
            TMPatlasName   = repmat(atlasName(a),1,Nthresh);
            end
            atlasFile      = cat(2, atlasFile(1:(a-1)),  TMPatlasFile,  atlasFile((a+1):end)  );
            isnorm         = cat(2, isnorm(1:(a-1)),     TMPisnorm,     isnorm((a+1):end)     );
            normthresh     = cat(2, normthresh(1:(a-1)), TMPnormthresh, normthresh((a+1):end) );
            issmry         = cat(2, issmry(1:(a-1)),     TMPissmry,     issmry((a+1):end)     );
            atlasName      = cat(2, atlasName(1:(a-1)),  TMPatlasName,  atlasName((a+1):end)  );
        else
            normthresh(a)  = normthreshSet{a};
            if ~dorename && isnorm(a)
            atlasName{a}   = sprintf('%s_%d',atlasName{a},round(normthresh(a)));
            end
        end
    end
else
% For non cell input
    atlasFile  = atlasFile{1};
    normthresh = normthreshSet{1};
    if ~dorename && isnorm
        atlasName   = sprintf('%s_%d',atlasName{1},round(normthresh(1)));
    else
        atlasName  = atlasName{1};
    end
end

        