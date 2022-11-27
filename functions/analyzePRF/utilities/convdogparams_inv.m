function    out = convdogparams_inv(x)

% out = convdogparams_inv(x)
%   convert parameters of DoG into [R C S G N SR GR] from [R C S1 G1 N S2 G2]
% 
% Arguments:
%   x = [R C S1 G1 N S2 G2] 	: parameters for DoG-CSS model
% 
% DoG = gaussian([R C S1=S*SS G1=G*GS N]) - gaussian([R C S2=S*SR*SS G2=G*GR*GS N]))
%   SR > 1      : sigma ratio for negative gaussian
%   GR = (0,1) 	: gain ratio for negative gaussian
%   SS          : sigma scale to convert one gaussian to center gaussian in DoG
%   GS          : gain scale to convert one gaussian to center gaussian in DoG
% 
% 
% See also: makegaussian2d, analyzePRFdog, modeldogcss, makedoggaussian, convdogparams

% 20221123 yuasa

%-- move the dimension of the parameters from 2 to 1
if iscolumn(x), x = x'; end
x = permute(x,[2,1,3:ndims(x)]);
sizx = size(x);
nx   = prod(sizx(2:end));
out  = zeros([7,sizx(2:end)]);

%-- convert for each DoG
for ii=1:nx
xi = x(:,ii);

GR = -xi(7)/xi(4);
SR = xi(6)/xi(3);

GS = 1/(1-GR);

%-- estimate initial SS with considering S2*SS is enough large
SS = sqrt(-1./(log(1+GR)./log(2)-1));
%-- estimate SS in loop to satisfy the equation about FWHM
dSS = 1; iter = 0;
while (abs(dSS)>eps && iter < 200)
    SS0=SS;
    SS = sqrt(-1./(log(1-GR*(1-2^(1-1/(SR*SS)^2)))./log(2)-1));
    dSS = (SS0-SS)/SS;
    iter = iter + 1;
end


out(:,ii)   = [xi(1), xi(2), xi(3)./SS, xi(4)./GS, xi(5), SR, GR];
end

%-- set the dimension of the parameters at 2
out = permute(out,[2,1,3:ndims(out)]);
