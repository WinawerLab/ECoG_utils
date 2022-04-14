function [I,varargout] = cellstrfind(C,patterns,expatterns,fullout)
% I      = CELLSTRFIND(cellstr,pattern)
% [I,I2] = CELLSTRFIND(cellstr,pattern)
% [I,I2] = CELLSTRFIND(cellstr,pattern,fulloutput)
%
% Returns index vector I at which cellstr contains the specified pattern,
% Pattern can be specified as a character vector or cellstr-array.
% Wildcard can be used in the pattern.
% I2 returns matching pairs between the specified cellstr and pattern.
% If fulloutput is true, I returns a logical vector instead of the index.
% 
% I      = CELLSTRFIND(cellstr,pattern,exclude_pattern)
% 
% You can exclude patterns specified in the exclude_pattern.
% If only the exclude_pattern is specified and the pattern is emptry, 
% it returns all index but supecified in the exclude_pattern.

% 160721: Yuasa
% 160811: Yuasa: enable multiple patterns
% 160906: Yuasa: bug fix (consider metacaracters as letters)
% 161219: Yuasa: Speed-up, add I2 output
% 190508: Yuasa: minor update
% 190729: Yuasa: add fulloutput option
% 190730: Yuasa: add exclude_pattern option
%                modify wildcard interpretation

narginchk(2,4);
nargoutchk(0,2);
if nargin<4||isempty(fullout),  fullout = false; end
if nargin<3||isempty(expatterns), expatterns = '';   end
if islogical(expatterns)||isnumeric(expatterns) % for (cellstr,pattern,fulloutput)
    if isscalar(expatterns), fullout = logical(expatterns); expatterns = '';  end
end
if ~isempty(expatterns)&&isempty(patterns), patterns = '*'; end % for only exclude
if isempty(C)||isempty(patterns), I = []; varargout{1} = []; return; end

% Check input 'C'
if iscellstr(C) || (exist('strings','builtin') && isstring(C))
    % string-array is also allowed ('exist' checks matlab ver)
    C = reshape(cellstr(C),1,[]);
else
    help cellstrfind
    error('First argument must be a cellstr.')
end

% Check input 'patterns' & resolve wildcard
try
    patterns = reshape(cellstr(patterns),1,[]);
    patterns = regexptranslate('wildcard', patterns);
catch
    error('''pettern'' and ''exclude_pattern'' must be a a character vector or cellstr-array.');
end

% find pattern
I_l  = zeros([numel(C), numel(patterns)]);
for l = 1:numel(patterns)
    pattern = ['^' patterns{l} '$'];
    I_l(:,l) = ~cellfun(@isempty,regexp(C,pattern,'once'));
end

% exclude
if ~isempty(expatterns)
    I_ex = cellstrfind(C,expatterns,true);
    I_l(I_ex,:) = false;
end

% output
if fullout, I  = logical(sum(I_l,2));
else,       I  = find(sum(I_l,2));
end
if fullout, varargout{1} = I_l;
else,       varargout{1} = I_l(I,:);
end

return
