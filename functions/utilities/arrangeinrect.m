function [nx,ny,nrate] = arrangeinrect(num,varargin)

% ARRANGEINRECT gives x and y length to arrange N items in a rectangle with
% the least blank panels. The priority is 1. to sutisfy the bouds of allow_rate, 
% 2. to minimize the number of blank panels, 3. to have the closest x/y rate 
% to the best_rate, 4. to have the closest x/y rate to the allow_rate(1).
% If it is impossible to arrange N with allow_rate, ARRANGEINRECT returns
% alternative pairs of x,y out of allow_rate, or reruns error for 'strict'
% mode. Here, the priority is 1. to have the closest x/y rate to the best_rate, 
% 2. to have the closest x/y rate to the allow_rate(1).
% 
% Usage: 
%   [x,y]     = arrangeinrect(N,best_rate,allow_rate)
%   [x,y,x/y] = arrangeinrect(N,best_rate,allow_rate)
%   [x,y,x/y] = arrangeinrect(N,best_rate,allow_rate,'strict')
%   [x,y,x/y] = arrangeinrect(N,best_rate,allow_rate,'silent')
% 
% Inputs: 
%   N:          numeric-array, numbers of items to arrange
%   best_rate:  number, desired rate of x/y (default = 1.0)
%   allow_rate: vector, bounds of x/y rate to find appropriate pairs of [x,y]
%                       (default = best_rate.*[1.0 2.0])
%               if allow_rate is scalar, the bounds set as [best_rate allow_rate]
% 

% 20170330 Yuasa
% 20190423 Yuasa: change usage of allow_rate
% 20200218 Yuasa: update for latest MATLAB version (numeric arrays are inhibited in warning)

% Dependency: error_backtraceoff

%-- set parameters
assert(nargin>=1,message('MATLAB:minrhs'));
assert(nargin<=5,message('MATLAB:maxrhs'));
assert(isnumeric(num),'MATLAB:arrangeinrect:nonNumericInput','''N'' must be a number');
% bestrate,allowrate,strict
isstrict = false;
issilent = false;
if nargin>1 && ischar(varargin{end}) 
    if strcmp(varargin{end},'strict')
        isstrict = true;
    elseif strcmp(varargin{end},'silent')
        issilent = true;
    end
end
if nargin>2 && ischar(varargin{end-1}) 
    if strcmp(varargin{end-1},'strict')
        isstrict = true;
    elseif strcmp(varargin{end-1},'silent')
        issilent = true;
    end
end
if (nargin - double(isstrict) - double(issilent))<2 || isempty(varargin{1})
    bestrate = 1.0;
else
    bestrate = varargin{1};
end
if isnumeric(bestrate)
    bestrate = bestrate(1);
else
    error('MATLAB:arrangeinrect:nonNumericInput','''best_rate'' must be a number')
end
if (nargin - double(isstrict) - double(issilent))<3 || isempty(varargin{2})
    allowrate = bestrate.*[1.0 2.0];
else
    allowrate = varargin{2};
end
if isnumeric(allowrate)
    if isscalar(allowrate)
        allowrate = [bestrate allowrate];
    else
        allowrate = reshape(allowrate(1:2),1,[]);
    end
else
    error('MATLAB:arrangeinrect:nonNumericInput','''allow_rate'' must be a number')
end

nx    = zeros(numel(num),1);
ny    = nx;     nrate    = nx;

dirc  = 2*(diff(allowrate)>=0)-1;

try
if issilent,    curwarn = warning('off');
else,           curwarn = warning('off','backtrace');
end

%%% warning for bad allow_rate
if bestrate<min(allowrate)||bestrate>max(allowrate)
    if isstrict
      error_backtraceoff(sprintf('best_rate (%g) is out of allow_rate ([%g %g])',bestrate,allowrate(1),allowrate(2)));
    else
      warning('best_rate (%g) is out of allow_rate ([%g %g])',bestrate,allowrate(1),allowrate(2));
    end
end
%%% main
for ilp=1:numel(num)

    begnx   = max(1,round(sqrt(num(ilp)*allowrate(1)))-dirc);
    endnx   = max(1,round(sqrt(num(ilp)*allowrate(2)))+dirc);

    %-- get mod list
    listnx   = reshape(begnx:dirc:endnx,[],1);
    listny   = ceil(num(ilp)./listnx);
      [listny,listidx]   = unique(listny,'stable');
      if dirc == -1,    listidx  = [reshape(listidx(2:end)-1,1,[]) length(listnx)];   end
      listnx    = listnx(listidx);
    listrate = listnx./listny;
    valididx = listrate>=min(allowrate)&listrate<=max(allowrate);
    if ~isempty(find(valididx, 1))
        listnx(~valididx)   = [];
        listny(~valididx)   = [];
        listrate(~valididx) = [];
        listmod  = mod(num(ilp),listny);
        
        %-- pick up the least-blank arrangement
        valididx = listmod==min(listmod);
        listnx(~valididx)   = [];
        listny(~valididx)   = [];
        listrate(~valididx) = [];
    elseif isstrict
        error_backtraceoff(sprintf('allow_rate ([%g %g]) is too strict for %g (i=%d)',allowrate(1),allowrate(2),num(ilp),ilp));
    else
        warning('allow_rate ([%g %g]) is too strict for %g (i=%d)',allowrate(1),allowrate(2),num(ilp),ilp);
    end
    %-- get nearest x/y rate pair
    diffrate = abs(listrate-bestrate);
    bestidx = find(diffrate == min(diffrate),1);
    %-- get output values
    nx(ilp)     = listnx(bestidx);
    ny(ilp)     = listny(bestidx);
    nrate(ilp)  = listrate(bestidx);
end

if nargout < 2
    nx = [nx ny];
end

warning(curwarn);
catch ME
    warning(curwarn);
    rethrow(ME);
end
