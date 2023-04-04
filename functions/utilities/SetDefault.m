function   varargout = SetDefault(varname, defval, varargin)
% SETDEFAULT assigns a value to the specified variable if the variable is
% not exist or empty.
% 
% Usage: 
%   issetdefault = SETDEFAULT(variable_name, default_value ,[keepempty, workspace])
%   issetdefault = SETDEFAULT(struct_name.field_name, default_value ,[keepempty, workspace])
% 
%   issetdefault = SETDEFAULT(variable_name, default_value ,[keepempty, workspace],'cell')
% 
% Inputs:
%   variable_name:  string, variable name.
%   (struct_name.field_name: strings, full set of structure name and field name)
%   default_value:  anything you want to assin to the variable or the structure filed.
%   keepempty:      false (default) or true, if true, SETDEFAULT assing a
%                   default value only when the variable is not exist.
%   workspace:      'base' or 'caller (default)', workspace in which to
%                   evaluate expression. See also EVALIN, ASSIGNIN.
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
% 20200220 Yuasa: add 'cell' mode
% 20220224 Yuasa: add 'base' mode
% 20230329 Yuasa: consider string as a kind of cell

narginchk(2,5)
issetdefault = false;

%-- arguments
if nargin < 3
    cellmode   = false;
    whichspace = 'caller';
    keepempty  = false;
else
    %-- keep empty
    charinputs = cellfun(@ischar,varargin);
    if sum(~charinputs) == 0
        keepempty = false;
    elseif sum(~charinputs) == 1
        keepempty = varargin{~charinputs};
    else
        noncharinputs = find(~charinputs);
        error('%s argument is invalid.', iptnum2ordinal(noncharinputs(2)+2));
    end
    %--  options
    opts = varargin(charinputs);
    valopts = ismember(opts,{'base','caller','cell'});
    assert(all(valopts),'''%s'' is invalid argument.\n',opts{~valopts});
    cellmode = any(ismember(opts,'cell'));
    if any(ismember(opts,'base')) && any(ismember(opts,'caller'))
        error('Some arguments are conflicted.');
    elseif any(ismember(opts,'base'))
        whichspace = 'base';
    else
        whichspace = 'caller';
    end
end

%-- keep empty
if keepempty
    varempty = false;
end

%-- cell mode
if cellmode && ~iscell(defval)
    defval = {defval};
end

%-- check existence
assert(ischar(varname),'1st input value must be a variable name.');
chkstrct = strsplit(varname,'.');
varexist = evalin(whichspace,sprintf('exist(''%s'',''var'')',chkstrct{1}));   % check existence
%-- check existence of subfield
if varexist
    if length(chkstrct) == 1    % need not to consider fields
        fldexst = true;
    else
        if evalin(whichspace,sprintf('isstruct(%s)||isempty(%s)',chkstrct{1},chkstrct{1}))    % check the structure
            getstrct = evalin(whichspace,sprintf('%s',chkstrct{1}));          % copy the structure
            [fldexst,nsubfld]  = issubfield(getstrct,chkstrct{2:end});
            if nsubfld==length(chkstrct) || ...
                    evalin(whichspace,sprintf('isstruct(%s)||isempty(%s)',...
                    strjoin(chkstrct(1:nsubfld),'.'),strjoin(chkstrct(1:nsubfld),'.')))
                eval(sprintf('%s = %s;',strjoin([{'getstrct'}, chkstrct(2:end)],'.'),'defval'));
                defval   = getstrct;
            else %-- skip if a part of fields exists but is not structure array
                fldexst  = true;
                varempty = false;
                cellmode  = false;
                varname  = strjoin(chkstrct(1:nsubfld),'.');
                wst = warning('backtrace','off');
                warning('''%s'' is not structure.',varname);
                warning(wst);
            end
        else
            fldexst  = true;
            varempty = false;
            cellmode  = false;
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
%     varempty = evalin(whichspace,sprintf('isempty(%s)',varname));             % check empty
% end
if varexist && fldexst
    SetDefault('varempty', evalin(whichspace,sprintf('isempty(%s)',varname)),'caller');             % check empty
end
if ~varexist || ~fldexst || varempty
    assignin(whichspace,chkstrct{1},defval);
    issetdefault = true;
elseif cellmode
    if evalin(whichspace,sprintf('~iscell(%s)&&~isstring(%s)',varname,varname))
        evalin(whichspace,sprintf('%s = {%s};',varname,varname));
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
