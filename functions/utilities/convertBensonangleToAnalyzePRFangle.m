function [channels] = convertBensonangleToAnalyzePRFangle(channels)
% converts the angle estimates from the benson14 atlas to the analyzePRF
% format.
%
% analyze PRF: 0 corresponds to the right horizontal meridian, 90 corresponds to the upper vertical meridian, and so on.
% benson14_angle: 0 corresponds to upper vertical meridian, 180 to lower vertical meridian
%
% [channels] = convertBensonangleToAnalyzePRFangle(channels);
%
% IG, 2020

% channels.benson14_angle = channels.benson14_angle+90; % THIS IS WRONG should take into account hemisphere!!!!

% From KEN:
%islefths = contains(channels.hemisphere, 'L');
%tmp = ((-1).^islefths).*channels.benson14_angle+90; % this gives negative angles??

if ~isfield(summary(channels), 'benson14_angle')
    error('no benson angle column present in this channel table')
end

if ~isfield(summary(channels), 'hemisphere')
    error('no hemisphere column present in this channel table')
end

% analyze PRF: 0 corresponds to the right horizontal meridian, 90 corresponds to the upper vertical meridian, and so on.
% benson14_angle: 0 corresponds to upper vertical meridian, 180 to lower vertical meridian
upper_right = contains(channels.hemisphere, 'L') & channels.benson14_angle <=90;
lower_right = contains(channels.hemisphere, 'L') & channels.benson14_angle >90;
left = contains(channels.hemisphere, 'R');
assert(sum(sum([upper_right lower_right left])) == height(channels))

tmp = channels.benson14_angle;

tmp(upper_right) = 90 - tmp(upper_right); % set right VF to 0-90, counterclockwise
tmp(lower_right) = (180 - tmp(lower_right)) + 240; % set right VF to 240-260, counterclockwise
tmp(left)        = tmp(left) + 90; % left VF to 90-240, counterclockwise

figure;hold on
subplot(1,2,1);
scatter(channels.benson14_angle(left), tmp(left), 100,'filled');
title('benson14_angle LH / RVF'); axis square
xlabel('original'); ylabel('converted');
subplot(1,2,2);
scatter(channels.benson14_angle(~left), tmp(~left), 100,'filled');
title('benson14_angle RH / LVF'); axis square
xlabel('original'); ylabel('converted');


channels.benson14_angle = tmp;

end