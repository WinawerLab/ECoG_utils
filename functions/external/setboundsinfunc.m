function varargout = setboundsinfunc(lbin,ubin,lbout,ubout,func,varargin)

% [output1,output2,...] = setboundsinfunc(lb_input,up_input,lb_output,ub_output,func,input1,input2,...)
%   
%   setboundsinfunc constrains inputs and outputs to satisfy the bounds.
%   
%   lb_input: Nx1 cell-array of lower boudns for N inputs
%   ub_input: Nx1 cell-array of upper boudns for N inputs
%   lb_output: Mx1 cell-array of lower boudns for M outputs
%   ub_output: Mx1 cell-array of upper boudns for M outputs
% 
% 
% [output1,output2,...] = setboundsinfunc(lb_input,up_input,lb_output,ub_output,func,input1,input2,...,Name,Value)
%   Following option is available
%   - 'OutputMode'  = 'constrain'(default), 'nan','zero'
%           'constrain': setboundsinfunc constrains inputs and outputs to satisfy the bounds
%           'nan':       setboundsinfunc outputs NaN when inputs or outputs exceed the bounds
%           'inf':       setboundsinfunc outputs Inf when inputs or outputs exceed the bounds
%           'zero':      setboundsinfunc outputs 0 when inputs or outputs exceed the bounds
%   

% 20200127 Yuasa
% 20200130 Yuasa: Upgrade <OutputMode>
% 
% Dependency: cellfind

%-- check options
narginchk(6,Inf);
numout    = nargout(func);
if numout>=0,   nargoutchk(0,numout);
else,           numout = max(abs(numout),nargout);   % case for using vargout in func
end
varargout = cell(1,numout);

%%%% OutputMode
optidx = cellfind(varargin,'OutputMode');
assert(isempty(optidx)||max(optidx)<length(varargin),'Invalid Option');
if ~isempty(optidx)
    OutputMode = lower(varargin{optidx(1)+1});
    assert(ischar(OutputMode)&&ismember(OutputMode,{'constrain','nan','inf','-inf','zero'}),'Invalid Option');
end
varargin([optidx, optidx+1]) = [];
    
%-- constrain inputs
outerflg_in = false;
numin = length(varargin);
if ~iscell(lbin),  lbin  = {lbin};   end
if ~iscell(ubin),  ubin  = {ubin};   end
lbin = [reshape(lbin,1,[]),cell(1,numin-numel(lbin))];
ubin = [reshape(ubin,1,[]),cell(1,numin-numel(ubin))];
%%% loop for inputs
for iin = 1:numin
    assert(isempty(lbin{iin})||isequal(size(lbin{iin}),size(varargin{iin})),'Bounds must be the same size as inputs');
    assert(isempty(ubin{iin})||isequal(size(ubin{iin}),size(varargin{iin})),'Bounds must be the same size as inputs');
    if ~isempty(lbin{iin}) && any(varargin{iin}<lbin{iin},'all')
        outerflg_in = true;
        varargin{iin}(varargin{iin}<lbin{iin}) = lbin{iin}(varargin{iin}<lbin{iin});
    end
    if ~isempty(ubin{iin}) && any(varargin{iin}>ubin{iin},'all')
        outerflg_in = true;
        varargin{iin}(varargin{iin}>ubin{iin}) = ubin{iin}(varargin{iin}>ubin{iin});
    end
end

%-- constrain outpus
if ~iscell(lbout), lbout = {lbout};   end
if ~iscell(ubout), ubout = {ubout};   end
lbout = [reshape(lbout,1,[]),cell(1,numout-numel(lbout))];
ubout = [reshape(ubout,1,[]),cell(1,numout-numel(ubout))];
%%% call function
[varargout{:}] = func(varargin{:});
%%% loop for outputs
for iout = 1:numout
    outerflg_out = false;
    assert(isempty(lbout{iout})||isequal(size(lbout{iout}),size(varargout{iout})),'Bounds must be the same size as outputs');
    assert(isempty(ubout{iout})||isequal(size(ubout{iout}),size(varargout{iout})),'Bounds must be the same size as outputs');
    if ~isempty(lbout{iout}) && any(varargout{iout}<lbout{iout},'all')
        outerflg_out = true;
        varargout{iout}(varargout{iout}<lbout{iout}) = lbout{iout}(varargout{iout}<lbout{iout});
    end
    if ~isempty(ubout{iout}) && any(varargout{iout}>ubout{iout},'all')
        outerflg_out = true;
        varargout{iout}(varargout{iout}>ubout{iout}) = ubout{iout}(varargout{iout}>ubout{iout});
    end
    if outerflg_out || outerflg_in
      switch OutputMode
        case 'nan',    varargout{iout} = nan(size(varargout{iout}));
        case 'inf',    varargout{iout} = inf(size(varargout{iout}));
        case '-inf',   varargout{iout} = -inf(size(varargout{iout}));
        case 'zero',   varargout{iout} = zeros(size(varargout{iout}));
      end
    end
end

