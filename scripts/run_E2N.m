% close all;
if isempty(which('electrode_to_nearest_node'))
    tbUse('ECoG_utils');
end

%% run script for a specific patient (with plotmesh on) to look at electrodes on brain

specs = [];

specs.pID           = 'som726'; % patient ID 
specs.atlasNames    =  {'benson14_varea', 'wang15_mplbl'};%{'glasser16_atlas'};%{'wang2015_atlas','wang15_mplbl','benson14_varea', 'benson14_eccen'}; 
                        % default is all of these maps: 
                        % {'wang2015_atlas', ...
                        % 'benson14_varea', ...
                        % 'benson14_eccen', ...
                        % 'benson14_angle', ...
                        % 'benson14_sigma', ...
                        % 'template_areas'};
                        % NOTE: including benson14_varea is required to
                        % be able to obtain benson14_eccen, angle and sigma
specs.plotmesh      = 'right';
specs.plotelecs     = 'yes';
specs.plotlabel     = 'yes';
%specs.fsDir         = '/Volumes/server/Projects/BAIR/Data/BIDS/visual_old/derivatives/freesurfer/';
out = electrode_to_nearest_node(specs,1);

%% to change one of the output figs to an 'empty' brain:

% Pick one of the figures and get its handle
fH = gcf;

% to change view to lateral:
%view(-90,0);
view(0,0);

% Find the patch object
h = findobj(fH, 'type', 'patch');

for ii = 1:length(h)
    % The background color is probably the modal color
    md = mode(h(ii).FaceVertexCData);
    
    % Set all vertices to the modal color
    h(ii).FaceVertexCData(:) = md;
    
    colorbar off
end

%% to generate a circle with a specific radius in Euclidean distance

radius = 25;

ii = 2; % do separately for each hemisphere (1 = left, 2 = right)

% If desired, get user to click a location
datacursormode on
dcmObj = datacursormode(fH);
set(dcmObj,'SnapToDataVertex','off','Enable','on', 'DisplayStyle', 'window')

waitforbuttonpress();
point1 = getCursorInfo(dcmObj);
x = point1.Position(1);
y = point1.Position(2);
z = point1.Position(3);


v = h(ii).Vertices;
xyzdiff = bsxfun(@minus, v, [x y z]);

dist = sqrt(sum(xyzdiff.^2,2));

mask = dist < radius;

h(ii).FaceVertexCData(mask) = 1;

%% views used on 10-15-18

view(0,0) %posterior
view(-90,0) % lateral LH
view(-30,0) % lateral LH rotated
% set(gcf, 'color','black')
