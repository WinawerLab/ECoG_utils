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
upper_RVF = contains(channels.hemisphere, 'L') & channels.benson14_angle <=90;
lower_RVF = contains(channels.hemisphere, 'L') & channels.benson14_angle >90;
leftVF = contains(channels.hemisphere, 'R');

assert(sum(sum([upper_RVF lower_RVF leftVF])) == height(channels)); % make sure we didn't miss any electrodes

converted = channels.benson14_angle;

converted(upper_RVF) = 90 - channels.benson14_angle(upper_RVF); % set right VF to 0-90, counterclockwise
converted(lower_RVF) = (180 - channels.benson14_angle(lower_RVF)) + 270; % set right VF to 240-360, counterclockwise
converted(leftVF)    = channels.benson14_angle(leftVF) + 90; % left VF to 90-270, counterclockwise

figure;hold on
subplot(1,2,1);
scatter(channels.benson14_angle(leftVF), converted(leftVF), 100,'filled');
title('benson14_angle LVF (RH)'); xlim([0 180]); ylim([0 360]); axis square
xlabel('original'); ylabel('converted');
subplot(1,2,2);
scatter(channels.benson14_angle(~leftVF), converted(~leftVF), 100,'filled');
title('benson14_angle RVF (LH)'); xlim([0 180]); ylim([0 360]);axis square
xlabel('original'); ylabel('converted');

channels.benson14_angle = converted;

end