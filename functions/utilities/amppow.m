function out = amppow(inp,P)

% C = AMPPOW(A,B) denotes element-by-element powers on the amplitude of A.
% AMPPOW(A,B) is equivalant to 
%   sign(A) .* abs(A).^B

% 20210617 Yuasa

out = sign(inp) .* (abs(inp) .^ P);