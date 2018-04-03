
function [out] = electrode_to_nearest_node(pID, thresh, plotmesh, plotlabel)

% pID = patient ID (number)
% thresh = maximal distance electrode to nearest node, in mm
% plotmesh = 'yes' or 'no'

% SETTINGS
if nargin < 4 || isempty(plotlabel)
    plotlabel = 'no';
end
if nargin < 3 || isempty(plotmesh)
    plotmesh = 'no';
end
if nargin < 2 || isempty(thresh)
    thresh = inf;% 20;
end
    
% RUN
disp(['running patient number ' num2str(pID)]);

% check which directory patient is in main BAIR directory, if not find it in SoM
patientDir = '/Volumes/server/Projects/BAIR/ECoG/';
if ~isdir([patientDir num2str(pID)])
    patientDir = '/Volumes/server/Projects/BAIR/ECoG/SoM/';
    if ~isdir([patientDir num2str(pID)])
        disp('patient directory not found - exiting');
        out = [];
        return
    end
end 

% read electrode coordinate file from BAIR/ECOG directory
D = dir([patientDir num2str(pID) '/*coor_T1*.txt']);
if ~isempty(D)
    elec_file = D(1).name;
    %disp(['reading ' D.name]);
    fid = fopen([patientDir num2str(pID) '/' elec_file]); E = textscan(fid, '%s%f%f%f%s'); fclose(fid);
    elec_xyz = [E{2} E{3} E{4}]; 
    elec_labels = E{1};
    elec_types = unique(E{5});
    disp(['types of electrodes found: ' elec_types{:}]);
else
    disp('no coordinate file found - exiting');
    out = [];
    return
end

% read surface reconstructions from Freesurfer_subjects directory
surf_file_rh = ['/Volumes/server/Freesurfer_subjects/som' num2str(pID) '/surf/rh.pial'];
surf_file_lh = ['/Volumes/server/Freesurfer_subjects/som' num2str(pID) '/surf/lh.pial'];
if exist(surf_file_rh, 'file') && exist(surf_file_lh, 'file')
    [vertices_r, faces_r] = read_surf(surf_file_rh);
    [vertices_l, faces_l] = read_surf(surf_file_lh);
    vertices = [vertices_r;vertices_l];
else
    disp('no freesurfer surfaces found - exiting');
    out = [];
    return
end

% match the electrode xyz with the nodes in the surfaces; find nearest
[indices, bestSqDist] = nearpoints(elec_xyz', vertices'); % function from vistasoft

% ignore electrodes that are more than [thresh] away from any surface node
keep_idx = find(bestSqDist<thresh);
indices = indices(keep_idx);
elec_xyz = elec_xyz(keep_idx,:);

% node indices to atlases
atlasNames = {'wang2015_atlas','template_areas'};

for a = 1:length(atlasNames)
    
    % get atlas label names
    if a == 1 % Wang atlas
        wangfile = '/Volumes/server/Projects/Kastner2015Atlas/ProbAtlas_v4/ROIfiles_Labeling.txt';
        fid = fopen(wangfile); textscan(fid, '%s', 5);  W = textscan(fid, '%s%f%s%s'); fclose(fid);
        area_labels = W{4};
    elseif a == 2 % Noah template
        area_labels = {'V1', 'V2', 'V3'}; 
    end

    % get atlases for this subject
    atlas_file_rh = ['/Volumes/server/Freesurfer_subjects/som' num2str(pID) '/surf/rh.' atlasNames{a} '.mgz'];
    atlas_file_lh = ['/Volumes/server/Freesurfer_subjects/som' num2str(pID) '/surf/lh.' atlasNames{a} '.mgz'];
    if exist(atlas_file_rh, 'file') && exist(atlas_file_lh, 'file')
        [atlas_rh] = load_mgh(atlas_file_rh);
        [atlas_lh] = load_mgh(atlas_file_lh);
    else
        disp('no atlas annotations found - exiting');
        out = [];
        return
    end
    
    % match nearest nodes to atlas labels
    atlas = [squeeze(atlas_rh);squeeze(atlas_lh)]; % concatenate hemis   
    atlas_elec = atlas(indices);
    [elec_indices] = find(atlas_elec);
    elec_labels_found = elec_labels(elec_indices);
    area_labels_found = area_labels(atlas_elec(elec_indices));
     
    % print to window how many maps were found, and
    % make a count per area (across hemispheres)
    area_count = zeros(1,length(area_labels));
    if ~isempty(elec_labels_found)
        disp(['found ' num2str(length(elec_labels_found)) ' electrodes in ' atlasNames{a} ': '])
        for i = 1:length(elec_labels_found)
            disp([elec_labels_found{i} ' in area ' area_labels_found{i}]);
            area_count(atlas_elec(elec_indices(i))) = area_count(atlas_elec(elec_indices(i))) + 1;
        end
    else 
        disp(['no electrodes in ' atlasNames{a}]);
    end
        
    % plot
    switch plotmesh
        case 'yes'
            figure('Name', num2str(pID)); hold on;
            
            plot_electrodes(elec_xyz, [1 1 1]*0.2,2);
            plot_electrodes(elec_xyz(elec_indices,:), [1 0 0],2);

            t_r = trimesh(faces_r+1, vertices_r(:,1), vertices_r(:,2), vertices_r(:,3), atlas_rh, 'FaceColor', 'flat'); 
            t_r.LineStyle = 'none';
            axis equal; hold on;
            t_l = trimesh(faces_l+1, vertices_l(:,1), vertices_l(:,2), vertices_l(:,3), atlas_lh, 'FaceColor', 'flat'); 
            t_l.LineStyle = 'none';
            cmap = [[1 1 1]*.7; colormap(flipud(jet(length(area_labels))))];
            colormap(cmap);
            caxis([0 max(atlas_rh)]);
            
            plot_electrodes(vertices(indices,:), [1 1 1]*0.8, 1);
            plot_electrodes(vertices(indices(elec_indices),:), [1 1 0], 1);
            switch plotlabel
                case 'yes'
                    for i = 1:size(elec_xyz,1)
                        [x, y, z] = adjust_elec_label(elec_xyz(i,:),2);
                        text('Position',[x y z],'String',elec_labels(i,:),'Color','w','VerticalAlignment','top');
                    end
            end
            axis off; set(gcf, 'color','black','InvertHardCopy', 'off');
            view(0,0);
            material dull;

            h=light; lightangle(h,  45, 45); lighting gouraud;
            h=light; lightangle(h, -45, 45); lighting gouraud;
            h=light; lightangle(h, -45, -90); lighting gouraud;
        otherwise
            % do not plot
    end
    
    % save 
    % also include count per atlas map to facilitate computing distribution
    % across patients?
    out.patientID = pID;
    out.(atlasNames{a}).elec_labels = elec_labels_found;
    out.(atlasNames{a}).area_labels = area_labels_found;
    out.(atlasNames{a}).area_count = area_count;
	out.(atlasNames{a}).area_names = area_labels;
    
end


end

function plot_electrodes(xyz, color, radius)
    [x, y, z] = sphere;
    if ~exist('radius', 'var'), radius = 2; end
        x = radius * x;
        y = radius * y;
        z = radius * z;
    for ii = 1:size(xyz,1)
       s = surface(xyz(ii,1)+x, xyz(ii,2)+y, xyz(ii,3)+z, 'FaceColor', color, 'EdgeColor', 'none', 'AmbientStrength',0.7);
       s.CDataMapping = 'direct';
    end
    %shading interp;
    %lighting gouraud;
    %material dull;
    %light;
end

function [x, y, z] = adjust_elec_label(xyz,radius)

    if ~exist('radius','var')
        radius = 2;
    end

    if xyz(1)>0
        x = xyz(1)+radius;
    else
        x = xyz(1)-radius;
    end

    if xyz(3)>0
        z = xyz(3)+radius;
    else
        z = xyz(3)-radius;
    end

    y = xyz(2);

end