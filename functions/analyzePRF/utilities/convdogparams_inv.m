function    out = convdogparams_inv(x)

% out = convdogparams_inv(x)
%   estimate equivalant parameters of [R C S1 G1 N S2 G2] from DoG parameters of [R C S G N SR GR] 
% 
% Arguments:
%   x = [R C S G N SR GR] 	: parameters for DoG-CSS model
% 
% DoG = gaussian([R C S1=S*SS G1=G*GS N]) - gaussian([R C S2=S*SR*SS G2=G*GR*GS N]))
%   SR > 1      : sigma ratio for negative gaussian
%   GR = (0,1) 	: gain ratio for negative gaussian
%   SS          : sigma scale to obtain the same width as the original positive gaussian
%   GS          : gain scale to obtain the same hight as the original positive gaussian
% 
% 
% See also: makegaussian2d, analyzePRFdog, modeldogcss, makedoggaussian

% hidden option
%   x = [R C S G N SR GR SS GS]

% 20191206 yuasa
% 20191213 yuasa: change preciseness for width
% under construction: currently work the same as convdogparams

%-- move the dimension of the parameters from 2 to 1
if iscolumn(x), x = x'; end
x = permute(x,[2,1,3:ndims(x)]);
sizx = size(x);
nx   = prod(sizx(2:end));
out  = zeros([7,sizx(2:end)]);

%-- convert for each DoG
for ii=1:nx
xi = x(:,ii);

if xi(6)==1 || xi(7)==0   % OG case
    xi(6) = 1;
    xi(7) = 0;
end
if numel(xi) < 9 || isnan(xi(9)),     GS = 1/(1-xi(7));
else,                               GS = xi(9);
end
if numel(xi) < 8 || isnan(xi(8))
    %-- estimate initial SS with considering S2*SS is enough large
    SS = sqrt(-1./(log(1+xi(7))./log(2)-1));
    %-- estimate SS in loop to satisfy the equation about FWHM
    dSS = 1; iter = 0;
    while (abs(dSS)>eps && iter < 200)
        SS0=SS;
        SS = sqrt(-1./(log(1-xi(7)*(1-2^(1-1/(xi(6)*SS)^2)))./log(2)-1));
        dSS = (SS0-SS)/SS;
        iter = iter + 1;
    end
end

%------ check DoG shape ------%
refS = (xi(3)./sqrt(posrect(xi(5))))*10;
xx   = linspace(-refS,refS,1e4);
yy   = posrect(xi(4)*GS).*(exp(-xx.^2./abs(xi(3)*SS).^2./2)./((2*pi*abs(xi(3)).^2))).^ posrect(xi(5)) ...
        - xi(7)*posrect(xi(4)*GS).*(exp(-xx.^2./abs(xi(6)*xi(3)*SS).^2./2)./((2*pi*abs(xi(3)).^2))).^ posrect(xi(5));
    
CxxIndex = yy > (max(yy)/3);
SxxIndex = [1:length(yy)] < find(yy==min(yy),1,'first') | [1:length(yy)] < find(yy==min(yy),1,'last');
%-----------------------------%

        
out(:,ii)   = [xi(1), xi(2), xi(3).*SS, xi(4).*GS, xi(5), xi(3).*SS.*xi(6), xi(4).*GS.*-xi(7)];
end

%-- set the dimension of the parameters at 2
out = permute(out,[2,1,3:ndims(out)]);
