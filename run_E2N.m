% close all;
if isempty(which('electrode_to_nearest_node'))
    tbUse('ECoG_utils');
end

%% run script for a specific patient (with plotmesh on) to look at electrodes on brain

specs = [];

specs.pID           = '692'; % patient ID number
specs.atlasNames    = {'wang2015_atlas', 'benson14_varea', 'benson14_eccen'}; 
                        % default is all of these maps: 
                        % {'wang2015_atlas', ...
                        % 'benson14_varea', ...
                        % 'benson14_eccen', ...
                        % 'benson14_angle', ...
                        % 'benson14_sigma', ...
                        % 'template_areas'};
                        % NOTE: including benson14_varea is required to
                        % be able to obtain benson14_eccen, angle and sigma
specs.plotmesh      = 'both'; % left, right, both, or none
specs.plotelecs     = 'no'; % yes or no

out = electrode_to_nearest_node(specs);

%% to change one of the output figs to an 'empty' brain:

% Pick one of the figures and get its handle
fH = gcf;

radius = 25;

% to change view to lateral:
view(-90,0);


% Find the patch object
h = findobj(fH, 'type', 'patch');

for ii = 1:length(h)
    % The background color is probably the modal color
    md = mode(h(ii).FaceVertexCData);
    
    % Set all vertices to the modal color
    h(ii).FaceVertexCData(:) = md;
    
    colorbar off
end

ii = 1;
% If desired, get user to click a location
datacursormode on
dcmObj = datacursormode(fH);
set(dcmObj,'SnapToDataVertex','off','Enable','on')

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



%% NEED TO FIX:UMCU

% UMCU patient electrode locations seem incorrect??

% specs = [];
% specs.pID      = 'beilen';
% specs.elecFile = '/Volumes/server/Projects/BAIR/Data/BIDS/visual/sub-umcubeilen/ses-umcuecogday03/ieeg/sub-beilen_ses-umcuecogday03_acq-clinical_electrodes.tsv';
% specs.fsDir    = '/Volumes/server/Freesurfer_subjects/umcu_beilen';
% specs.thresh   = []; 
% specs.patientPool = 'BAIR';
