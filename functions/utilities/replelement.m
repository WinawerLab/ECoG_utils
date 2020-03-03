function datamat = replelement(datamat, idx, vals)

% B = REPLELEMENT(A,idx,p)
%   replaces elements in matrix A indicated by idx with p.
% 
%   replelement(A,idx,p) is equivalant to A(idx) = p, where
%   p is a scalar or a matrix which has the same size as A(idx).

% 20200129 Yuasa

narginchk(3,3);

datamat(idx) = vals;