function TF = istablefield(T,fields)

% ISTABLEFIELD(T,FIELD) works as ISFIELD for TABLE.
% 
% Usage: 
%   ISTABLEFIELD(T,FIELD);
%   TF = ISTABLEFIELD(T,FIELDNAMES);
% 
% See also ISFIELD, TABLE.

% 20191105 Yuasa

assert(istable(T),'''%s'' is not table',inputname(1));
TF = ismember(fields, T.Properties.VariableNames);