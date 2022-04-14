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
% KY, 2022 - add error message

if ~isfield(summary(channels), 'benson14_angle')
    error('no benson angle column present in this channel table')
end

if ~isfield(summary(channels), 'hemisphere')
    error('no hemisphere column present in this channel table')
end

% analyze PRF: 0 corresponds to the right horizontal meridian, 90 corresponds to the upper vertical meridian, and so on.
% benson14_angle: 0 corresponds to upper vertical meridian, 180 to lower vertical meridian
rightVF = contains(channels.hemisphere, 'L');
leftVF  = contains(channels.hemisphere, 'R');

assert(all(or(rightVF,leftVF)),'Failed to specify left/right hemisphere in some electrodes'); % make sure we didn't miss any electrodes

% set right VF to 0-90, 270-360, and left VF to 90-270, counterclockwise
converted = mod( (-1).^(rightVF) .* channels.benson14_angle + 90, 360 );

% % debug
% figure;hold on
% subplot(1,2,1);
% scatter(channels.benson14_angle(leftVF), converted(leftVF), 100,'filled');
% title('benson14_angle LVF (RH)'); xlim([0 180]); ylim([0 360]); axis square
% xlabel('original'); ylabel('converted');
% subplot(1,2,2);
% scatter(channels.benson14_angle(~leftVF), converted(~leftVF), 100,'filled');
% title('benson14_angle RVF (LH)'); xlim([0 180]); ylim([0 360]);axis square
% xlabel('original'); ylabel('converted');

channels.benson14_angle = converted;

end