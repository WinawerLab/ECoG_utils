function I = cellfind(C,pattern)
% I = cellfind(cellArray,[pattern])
% I = cellfind(cellArray,cellArray)
%
% Returns index vector I of cells containing a desired pattern. Pattern can
% be anything that would be in a cell: a num, string, struct, another cell,
% etc. If pattern is omitted, I returns the indices of non-empty cells.
% 
% Multiple patterns can be specified as a cell-array. Then I returns the
% indices which match either of the patterns.
%
% Note the algorithm for finding a pattern isn't
% optimized to be fast, just convenient.

% ras 11/04: updated to include pattern finding
% bw/ab    : caught empty Cell condition
% 180510: Yuasa: modified

if isempty(C), I = []; return; end
% if ~ismatrix(C), error('Cellfind doesn''t work on arrays > 2 dimensions.'); end

if ~iscell(C)
    help cellfind
    error('First argument must be a cell.')
end

I=zeros(size(C));
if nargin < 2
    % no pattern, just find non-empty cells
	
	for t = 1:size(C,1)
        for u = 1:size(C,2)
            a = C{t,u};
            I(t,u)=~isempty(a);
        end
	end
else
    if ~iscell(pattern), pattern = {pattern};   end    
    % find pattern
    
    for iP = 1:numel(pattern)
        for iC = 1:numel(C)
            I(iC) = I(iC) + isequal(C{iC},pattern{iP});
        end
    end
end

I = find(I);

return
