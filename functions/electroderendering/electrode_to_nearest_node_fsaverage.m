
function [out] = electrode_to_nearest_node_fsaverage(specs)

% pID = patient ID (number)
% thresh = maximal distance electrode to nearest node, in mm
% plotmesh = 'yes' or 'no'

if ~isfield(specs, 'atlasNames') || isempty(specs.atlasNames)
    specs.atlasNames = {'wang2015_atlas','wang15_mplbl', 'benson14_varea', 'benson14_eccen', 'benson14_angle', 'benson14_sigma', 'template_areas'};
end

pID    = specs.pID;
thresh = specs.thresh;

if isempty(specs.thresh)
    thresh = inf;
end

if ~isfield(specs, 'fsDir') || isempty(specs.fsDir)
    specs.fsDir = fullfile(filesep, 'Volumes', 'server', 'Freesurfer_subjects'); 
end

plotlabel = specs.plotlabel;
plotmesh  = specs.plotmesh;

% RUN
disp(['running patient ' num2str(pID)]);

% Read electrode coordinate file from BAIR/ECOG directory
%elec_file = specs.elecFile;

switch specs.patientPool
    case 'BAIR'
        if exist(elec_file, 'file')

            % Prefer to use readtable for tsv files because it doesn't require
            % knowing the order of the columns beforehand (as textscan does). 
            disp(['reading ' elec_file]); % if there are multiple coor_T1 files for separate hemisphere, D(1) will always be the full list
            E = readtable(elec_file, 'FileType', 'text');
            %elec_labels = E.name;
            if iscell(E.x)
                for ii = 1:height(E)
                    elec_xyz(ii,1) = str2double(E.x{ii});
                    elec_xyz(ii,2) = str2double(E.y{ii});
                    elec_xyz(ii,3) = str2double(E.z{ii});    
                end
            else
                elec_xyz = [E.x E.y E.z];
            end
        else
            disp('electrode coordinate file not found - exiting. (NYU: is the server mapped?)');
            out = [];
            return
        end
    case 'SOM'
        
        % check which directory patient is in main BAIR directory, if not find it in SoM
        patientDir = '/Volumes/server/Projects/BAIR/Data/Raw/ECoG/';
        if ~isdir([patientDir num2str(pID)])
            patientDir = '/Volumes/server/Projects/BAIR/Data/Raw/ECoG/SoM/';
            if ~isdir([patientDir num2str(pID)])
                disp('patient directory not found - exiting');
                out = [];
                return
            end
        end
        
        % read electrode coordinate file from ECOG/SoM directory
        D = dir([patientDir num2str(pID) '/*coor_MNI*.txt']);
        if ~isempty(D)
            elec_file = D(1).name;
            disp(['reading ' D(1).name]); % if there are multiple coor_T1 files for separate hemisphere, D(1) will always be the full list
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
end

% Read surface reconstructions from Freesurfer_subjects directory
surf_file_rh = ['/Applications/freesurfer/subjects/fsaverage/surf/rh.pial'];
surf_file_lh = ['/Applications/freesurfer/subjects/fsaverage/surf/lh.pial'];

%surf_file_rh = fullfile(specs.fsDir, sprintf('som%s', specs.pID), 'surf', 'rh.pial');
%surf_file_lh = fullfile(specs.fsDir, sprintf('som%s', specs.pID), 'surf', 'lh.pial');

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
%atlasNames = {'wang2015_atlas'};%;,'template_areas', 'template_angle', 'template_eccen'};

for a = 1:length(specs.atlasNames)
    
    currentAtlas = specs.atlasNames{a};
    
    atlas_version = 'v4_0';
    
    % Get names / colormaps associated with each atlas
    switch currentAtlas

        case {'wang2015_atlas', 'wang15_mplbl'} % Wang atlas

            % Labels come from: '/Volumes/server/Projects/Kastner2015Atlas/ProbAtlas_v4/ROIfiles_Labeling.txt';
            area_labels =  {'V1v','V1d', ...
                            'V2v','V2d', ...
                            'V3v','V3d', ...
                            'hV4', ...
                            'VO1','VO2', ...
                            'PHC1','PHC2', ...
                            'TO2','TO1', ...
                            'LO2','LO1', ...
                            'V3b','V3a', ...
                            'IPS0','IPS1','IPS2','IPS3','IPS4','IPS5', ...
                            'SPL1','FEF'}; 
            area_cmap   = [102  51   0; 255   0   0; 
                           102   0  51; 255 128   0;
                            51   0 102; 255 255   0;
                             0 128 255;
                             0 102  51; 153 153   0;
                           204 102   0; 204   0   0;
                             0  76 153; 153  76   0;
                           255  51 153; 102   0 204;
                             0 255 255;   0 255   0;
                           255 153 153; 255 204 153; 255 255 153; 153 255 153; 153 255 255; 153 153 255; 
                           255 153 255; 255 178 102]./255;
            out.(currentAtlas).area_names   = area_labels;
            atlas_version = 'v1_0';
            
        case 'benson14_varea' % Noah template: areas 

            % Labels come from Noah email:  values in order from 1-12: V1 V2 V3 hV4 VO1 VO2 LO1 LO2 TO1 TO2 V3b V3a
            area_labels = {'V1', ...
                           'V2', ...
                           'V3', ...
                           'hV4', ...
                           'VO1', 'VO2', ...
                           'LO1', 'LO2', ...
                           'TO1', 'TO2', ...
                           'V3b', 'V3a'}; 
            area_cmap   = [255   0   0; 
                           255 128   0; 
                           255 255   0; 
                             0 128 255;
                             0  76 153; 153  76   0;
                             0 102  51; 153 153   0;  
                           255  51 153; 102   0 204;
                             0 255 255;   0 255   0]./255;
            out.(currentAtlas).area_names   = area_labels;

        case 'template_areas' % old Noah templates

            area_labels = {'V1', 'V2', 'V3'}; 
            area_cmap   = [255 255 0; 0 255 255; 0 0 255]./255;
            out.(currentAtlas).area_names   = area_labels;

        otherwise

            switch specs.plotmesh
                case {'both', 'left', 'right'}

                    fprintf('[%s] Loading colormap_%s \n',mfilename, currentAtlas);
                    A = load(['colormap_' currentAtlas]);

                    switch currentAtlas
                        case 'benson14_eccen'
                            area_cmap = A.cmap(15:end-200,2:4);
                        otherwise
                            area_cmap = A.cmap(:,2:4); 
                            % actually different colormaps should be
                            % used for left and right polar angle, but
                            % not possible to implement for 'both
                            % hemisphere ' plot
                    end

                    %atlas_range = range(atlas);
                    %cmap_index = round(linspace(1,length(area_cmap),atlas_range));
                    %area_cmap = area_cmap(cmap_index,:);   
            end
    end

    
    atlasDir = '/Users/winawerlab/Library/Python/2.7/lib/python/site-packages/neuropythy/lib/data/fsaverage/surf/';
    
%     % get atlases for this subject
    atlas_file_rh = fullfile(atlasDir,  sprintf('rh.%s.%s.mgz', currentAtlas, atlas_version));
	atlas_file_lh = fullfile(atlasDir,  sprintf('lh.%s.%s.mgz', currentAtlas, atlas_version));

    if exist(atlas_file_rh, 'file') && exist(atlas_file_lh, 'file')
        [atlas_rh] = load_mgh(atlas_file_rh);
        [atlas_lh] = load_mgh(atlas_file_lh);
    else
        fprintf('[%s] No annotations found for %s \n',mfilename, currentAtlas);
        out.(currentAtlas) = [];
        continue % move on to the next atlas
    end
    
%     % match nearest nodes to atlas labels
%     atlas = [squeeze(atlas_rh);squeeze(atlas_lh)]; % concatenate hemis   
%     atlas_elec = atlas(indices);
%     [elec_indices] = find(atlas_elec);
%     elec_labels_found = elec_labels(elec_indices);
%      
%     % print to window how many maps were found, and
%     % make a count per area (across hemispheres)
%     area_count = zeros(1,length(area_labels));
%     if a < 3
%         area_labels_found = area_labels(round(atlas_elec(elec_indices)));
%         if ~isempty(elec_labels_found)
%             disp(['found ' num2str(length(elec_labels_found)) ' electrodes in ' atlasNames{a} ': '])
%             for i = 1:length(elec_labels_found)
%                 disp([elec_labels_found{i} ' in area ' area_labels_found{i}]);
%                 area_count(atlas_elec(elec_indices(i))) = area_count(atlas_elec(elec_indices(i))) + 1;
%             end
%         else 
%             disp(['no electrodes in ' atlasNames{a}]);
%         end
%     end
%         
    % plot
    switch plotmesh
        case {'left', 'right', 'both'}
            switch pID 
                case '648'
                    
                    %figure(1)%,'Name', ['fsaverage ' atlasNames{a}]); hold on;
                    figure('Name', 'fsaverage '); hold on;
                    t_r = trimesh(faces_r+1, vertices_r(:,1), vertices_r(:,2), vertices_r(:,3), atlas_rh, 'FaceColor', 'flat');
                    t_r.LineStyle = 'none';
                    axis equal; hold on;
                    t_l = trimesh(faces_l+1, vertices_l(:,1), vertices_l(:,2), vertices_l(:,3), atlas_lh, 'FaceColor', 'flat');
                    t_l.LineStyle = 'none';
                    cmap = [[1 1 1]*.7; colormap(flipud(jet(length(area_labels))))];
                    %cmap = [1 1 1]*.7;
                    colormap(cmap);
                    caxis([0 max(atlas_rh)]);
                    
                    plot_electrodes(elec_xyz, [1 1 1]*0.2,2);
                    %plot_electrodes(elec_xyz(elec_indices,:), [0 0 0],2);
                    %plot_electrodes(vertices(indices,:), [1 1 1]*0.8, 1);
                    %plot_electrodes(vertices(indices(elec_indices),:), [1 1 1], 1);
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
                    plot_electrodes(elec_xyz, [1 1 1]*0.2,2);
                    %plot_electrodes(elec_xyz(elec_indices,:), [0 0 0],2);
                    %plot_electrodes(vertices(indices,:), [1 1 1]*0.8, 1);
                    %plot_electrodes(vertices(indices(elec_indices),:), [1 1 1], 1);
            end
            
       
        otherwise
            % do not plot
    end
    axis tight
    % save 
    % also include count per atlas map to facilitate computing distribution
    % across patients?
    out.patientID = pID;
    %out.(atlasNames{a}).elec_labels = elec_labels_found;
    %out.(atlasNames{a}).area_labels = area_labels_found;
    %out.(atlasNames{a}).area_count = area_count;
	%out.(atlasNames{a}).area_names = area_labels;
    
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