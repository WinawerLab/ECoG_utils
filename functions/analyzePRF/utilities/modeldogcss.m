function    out = modeldogcss(x,y,stimulus,res,varargin)

% out = modeldogcss(x,y,stimulus,res [,xx,yy,ang,omitexp,noneggain])
%   DoG CSS model responses for stimulus
% 
% Arguments:
%   x = [R C S G N]     : parameters for CSS model
%   y = [SR GR]         : parameters for DoG model
%   stimulus            : vector of stimlus image
% 
% DoG = gaussian([R C S*SS G*GS N]) - gaussian([R C S*SR*SS G*GR*GS N]))
%   SR > 1      : sigma ratio for negative gaussian
%   GR = (0,1) 	: gain ratio for negative gaussian
%   SS          : sigma scale to convert one gaussian to center gaussian in DoG
%   GS          : gain scale to convert one gaussian to center gaussian in DoG
% 
% Example:
%   modelfun = @(pp,dd) conv2run(makedoggaussian(pp(1:5),pp(6:end),dd,res,xx,yy,0,0),hrf,dd(:,prod(res)+1));
% 
% See also: makegaussian2d, analyzePRFdog, makedoggaussian, convdogparams

% hidden option
%   y = [SR GR SS GS nanoutput]
%   If nanoutput is true, then it outputs NaN image when SR or GR is out of threshold.
%   y = [SR GR] is equivalant to [SR GR NaN NaN false]

% Dependency: amppow

% 20191119 yuasa
% 20191126 yuasa: enable nanoutput as hidden option
% 20191212 yuasa: comment out limitation for sigma ratio
% 20191213 yuasa: change preciseness for width
% 20210616 yuasa: minor change to allow negative values for stimulus
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

if nargin<4 || isempty(res)
    res    = round(sqrt(length(stimulus)-1));
    if length(stimulus)~=res^2+1
        narginchk(4,inf);
    end
end 
resmx  = max(res);
outofbnds = false;
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
    out = nan(size(stimulus,1),1);
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
    % estimate initial SS with considering SR*SS is enough large
    SS = sqrt(-1./(log(1+y(2))./log(2)-1));
    % estimate SS in loop to satisfy the equation about FWHM
    dSS = 1; iter = 0;
    while (abs(dSS)>eps && iter < 200)
        SS0=SS;
        SS = sqrt(-1./(log(1-y(2)*(1-2^(1-1/(y(1)*SS)^2)))./log(2)-1));
        dSS = (SS0-SS)/SS;
        iter = iter + 1;
    end
else
    SS = y(3);
end

if noneggain,  gainfnc = @(x)posrect(x);
else,          gainfnc = @(x)x;
end
pp   = [[x(1), x(1)];
        [x(2), x(2)];
        x(3).*SS.*[1 y(1)];
        gainfnc(x(4)).*GS.*[1 -y(2)];
        posrect([x(5), x(5)])];
    
% if abs(pp(3,2))/sqrt(pp(5,2)) > resmx*1.5  % pRF size is enough larger than visual field
%     out = nan(size(stimulus,1),1);
%     return;
% end

out    = amppow((stimulus*...
            [[vflatten(placematrix(zeros(res),makegaussian2d(resmx,pp(1,1),pp(2,1),abs(pp(3,1)),abs(pp(3,1)),varargin{:}) / (2*pi*abs(x(3))^2)));0 ],...
             [vflatten(placematrix(zeros(res),makegaussian2d(resmx,pp(1,2),pp(2,2),abs(pp(3,2)),abs(pp(3,2)),varargin{:}) / (2*pi*abs(x(3))^2)));0 ]])...
            , posrect(pp(5,:))) * pp(4,:)';

%%% for debug
% xx = (-res:res)./2;
% figure; plot(xx,posrect(x(4))*exp(-xx.^2./(2*x(3).^2))./(2*pi*abs(x(3))^2));
% hold on; plot(xx,pp(4,1)*exp(-xx.^2./(2*pp(3,1).^2))./(2*pi*abs(x(3))^2)+pp(4,2)*exp(-xx.^2./(2*pp(3,2).^2))./(2*pi*abs(x(3))^2));
% figure; plot(xx,pp(4,1)*exp(-xx.^2./(2*pp(3,1).^2))./(2*pi*abs(x(3))^2),xx,pp(4,2)*exp(-xx.^2./(2*pp(3,2).^2))./(2*pi*abs(x(3))^2));
