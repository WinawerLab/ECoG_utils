function h = ciplot(lower,upper,x,colour, transp);
     
% ciplot(lower,upper)       
% ciplot(lower,upper,x)
% ciplot(lower,upper,x,colour)
% ciplot(lower,upper,x,colour,transp)
%
% Plots a shaded region on a graph between specified lower and upper confidence intervals (L and U).
% l and u must be vectors of the same length.
% Uses the 'fill' function, not 'area'. Therefore multiple shaded plots
% can be overlayed without a problem. Make them transparent for total visibility.
% x data can be specified, otherwise plots against index values.
% colour can be specified (eg 'k'). Defaults to blue.
% Raymond Reynolds 24/11/06
 
if length(lower)~=length(upper)
    error('lower and upper vectors must be same length')
end
 
if nargin < 5 
    transp = 0;
end
 
if nargin<4 || isempty(colour)
    colour='k';
end
 
if nargin<3 || isempty(x)
    x=1:length(lower);
end
 
% convert to row vectors so fliplr can work
if find(size(x)==(max(size(x))))<2
x=x'; end
if find(size(lower)==(max(size(lower))))<2
lower=lower'; end
if find(size(upper)==(max(size(upper))))<2
upper=upper'; end
 
h = fill([x fliplr(x)],[upper fliplr(lower)],colour);
set(h, 'EdgeColor', 'none');
if transp > 0 
    %camlight; lighting flat; 
    set(h,'facealpha',transp);
end


 