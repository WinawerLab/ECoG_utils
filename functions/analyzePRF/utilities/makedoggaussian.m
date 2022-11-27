function    out = makedoggaussian(x,y,res,varargin)

% out = makedoggaussian(x,y,res [,xx,yy,ang,omitexp])
%   makes 2D gaussian image based on DoG CSS model
% 
% Arguments:
%   x = [R C S G N]     : parameters for CSS model
%   y = [SR GR]         : parameters for DoG model
% 
% DoG = gaussian([R C S*SS G*GS N]) - gaussian([R C S*SR*SS G*GR*GS N]))
%   SR > 1      : sigma ratio for negative gaussian
%   GR = (0,1) 	: gain ratio for negative gaussian
%   SS          : sigma scale to convert one gaussian to center gaussian in DoG
%   GS          : gain scale to convert one gaussian to center gaussian in DoG
% 
% Example:
%   imagesc(makedoggaussian(pp(1:5),pp(6:end),res,xx,yy,0,0)));
% 
% See also: makegaussian2d, analyzePRFdog, modeldogcss, convdogparams

% hidden option
%   y = [SR GR SS GS nanoutput]
%   If nanoutput is true, then it outputs NaN image when SR or GR is out of threshold.
%   y = [SR GR] is equivalant to [SR GR NaN NaN false]

% Dependency: amppow

% 20191119 yuasa
% 20191126 yuasa: enable nanoutput as hidden option
% 20191213 yuasa: change preciseness for width
% 20221123 Yuasa: enable negative gain

if length(varargin)>4
    noneggain = varargin{5};
    varargin(5:end) = [];
else
    noneggain = [];
end
if isempty(noneggain)
    noneggain = false;
end

resmx  = max(res);
outofbnds = false;
% if y(1)>10
%     warning('sigma ratio should be smaller than 10');
%     outofbnds = true;
% end
if y(1)<1
    warning('sigma ratio must be lager than 1');
    outofbnds = true;
    y(1)=1;
end
if y(2)<0 || y(2)>1
    warning('gain ratio must be lager than 0 and smaller than 1');
    outofbnds = true;
    if y(2)<0,  y(2)=0;
    else,       y(2)=1;
    end
end

if numel(y) >= 5 && y(5) && outofbnds
    out = nan(res);
    return;
end
if y(1)==1 || y(2)==0   % OG case
    y(1) = 1;
    y(2) = 0;
end
if numel(y) < 4 || isnan(y(4)),     GS = 1/(1-y(2));
else,                               GS = y(4);
end
if numel(y) < 3 || isnan(y(3))
    % estimate initial SS with considering S2*SS is enough large
    SS = sqrt(-1./(log(1+y(2))./log(2)-1));
    % estimate SS in loop to satisfy the equation about FWHM
    dSS = 1; iter = 0;
    while (abs(dSS)>eps && iter < 200)
        SS0=SS;
        SS = sqrt(-1./(log(1-y(2)*(1-2^(1-1/(y(1)*SS)^2)))./log(2)-1));
        dSS = (SS0-SS)/SS;
        iter = iter + 1;
    end
end

if noneggain,  gainfnc = @(x)posrect(x);
else,          gainfnc = @(x)x;
end
pp   = [[x(1), x(1)];
        [x(2), x(2)];
        x(3).*SS.*[1 y(1)];
        gainfnc(x(4)).*GS.*[1 -y(2)];
        posrect([x(5), x(5)])];

gauss1 = pp(4,1) * amppow(placematrix(zeros(res),...
            makegaussian2d(resmx,pp(1,1),pp(2,1),abs(pp(3,1)),abs(pp(3,1)),varargin{:}) / (2*pi*abs(x(3))^2))...
            , pp(5,1));
gauss2 = pp(4,2) * amppow(placematrix(zeros(res),...
            makegaussian2d(resmx,pp(1,2),pp(2,2),abs(pp(3,2)),abs(pp(3,2)),varargin{:}) / (2*pi*abs(x(3))^2))...
            , pp(5,2));
        
out    = gauss1 + gauss2;
