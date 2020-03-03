function   varargout = SetDefault(varname, defval, keepempty, chkcell)
% SETDEFAULT assigns a value to the specified variable if the variable is
% not exist or empty.
% 
% Usage: 
%   issetdefault = SETDEFAULT(variable_name, default_value ,[keepempty])
%   issetdefault = SETDEFAULT(struct_name.field_name, default_value ,[keepempty])
% 
%   issetdefault = SETDEFAULT(variable_name, default_value ,[keepempty],'cell')
% 
% Inputs:
%   variable_name:  string, variable name.
%   (struct_name.field_name: strings, full set of structure name and field name)
%   default_value:  anything you want to assin to the variable or the structure filed.
%   keepempty:      false (default) or true, if true, SETDEFAULT assing a
%                   default value only when the variable is not exist.
% 
%   If 'cell' is specified, SETDEFAULT modifies the variable as cell-array
%   or assigns a value as cell-array.
% 
% Output:
%   issetdefault: true or false, if SETDEFAULT assigned a default value
% 
% Examples:
%   SetDefault('SampleRate',512);
%   SetDefault('cfg.fsample',512);
% 

% 20170622 Yuasa
% 20200220 Yuasa: cell mode

narginchk(2,4)

if nargin < 3
    chkcell = false;
    keepempty = false;
elseif nargin < 4 && ischar(keepempty)
    chkcell   = keepempty;
    keepempty = false;
else
    keepempty = logical(keepempty);
    if nargin < 4,    chkcell = false;    end
end
if ischar(chkcell)
    assert(strcmpi(chkcell,'cell'),'''%s'' is invalid input',chkcell);
    chkcell   = true;
else
    chkcell   = logical(chkcell);
end

issetdefault = false;

%-- keep empty
if keepempty
    varempty = false;
end

%-- cell mode
if chkcell && ~iscell(defval)
    defval = {defval};
end

%-- check existence
assert(ischar(varname),'1st input value must be a variable name.');
chkstrct = strsplit(varname,'.');
varexist = evalin('caller',sprintf('exist(''%s'',''var'')',chkstrct{1}));   % check existence
%-- check existence of subfield
if varexist
    if length(chkstrct) == 1    % need not to consider fields
        fldexst = true;
    else
        if evalin('caller',sprintf('isstruct(%s)||isempty(%s)',chkstrct{1},chkstrct{1}))    % check the structure
            getstrct = evalin('caller',sprintf('%s',chkstrct{1}));          % copy the structure
            [fldexst,nsubfld]  = issubfield(getstrct,chkstrct{2:end});
            if nsubfld==length(chkstrct) || ...
                    evalin('caller',sprintf('isstruct(%s)||isempty(%s)',...
                    strjoin(chkstrct(1:nsubfld),'.'),strjoin(chkstrct(1:nsubfld),'.')))
                eval(sprintf('%s = %s;',strjoin([{'getstrct'}, chkstrct(2:end)],'.'),'defval'));
                defval   = getstrct;
            else %-- skip if a part of fields exists but is not structure array
                fldexst  = true;
                varempty = false;
                chkcell  = false;
                varname  = strjoin(chkstrct(1:nsubfld),'.');
                wst = warning('backtrace','off');
                warning('''%s'' is not structure.',varname);
                warning(wst);
            end
        else
            fldexst  = true;
            varempty = false;
            chkcell  = false;
            varname  = chkstrct{1};
            wst = warning('backtrace','off');
            warning('''%s'' is not structure.',varname);
            warning(wst);
        end
    end
elseif length(chkstrct) ~= 1
    getstrct = [];
    eval(sprintf('%s = %s;',strjoin([{'getstrct'}, chkstrct(2:end)],'.'),'defval'));
    defval   = getstrct;
end
%-- check empty
% if ~exist('varempty','var') && varexist && fldexst
%     varempty = evalin('caller',sprintf('isempty(%s)',varname));             % check empty
% end
if varexist && fldexst
    SetDefault('varempty', evalin('caller',sprintf('isempty(%s)',varname)));             % check empty
end
if ~varexist || ~fldexst || varempty
    assignin('caller',chkstrct{1},defval);
    issetdefault = true;
elseif chkcell
    if evalin('caller',sprintf('~iscell(%s)',varname))
        evalin('caller',sprintf('%s = {%s};',varname,varname));
    end
end
if nargout~=0
    varargout = {issetdefault};
end


function [A,k] = issubfield(strct, varargin)
% ISSUBFIELD try ISFIELD recursively

if isempty(strct)||nargin<2
  A = false;
  k = 1;
else
  fieldname = varargin;
  A = true;
  for k = 1:numel(fieldname)
    if isfield(strct, fieldname{k})
      strct = strct.(fieldname{k});
    else
      A = false;
      return;
    end
  end
  k = k+1;
end
