
function [out] = electrode_to_nearest_node(specs)

% This function matches a list of electrode locations in an ECoG patient to
% the nearest node in their T1s freesurfer pial surface reconstruction,
% determines which of those nodes fall within a set of visual regions
% specified by two probabilistic atlases, and extracts node information.
% Note that electrode coordinates should be specified in the T1 volume
% space that was used to obtain the freesurfer reconstruction.
%
% INPUT: a struct containing the following fields:
% specs.pID         = patient ID 
% specs.elecFile    = file with electrode coordinates (include full path)
% specs.fsDir       = freesurfer directory of patient (include full path);
%                     should contain the wang and benson atlases that can
%                     be obtained through the nben/neuropythy docker
% specs.thresh      = maximum allowed distance between electrode and node,
%                     in mm (if left empty, thresh is infinite, meaning
%                     that the electrode can be infinitely far)
% specs.patientPool = should be either 'BAIR' or 'SOM' (may be removed
%                     later but necessary now to deal with differences in
%                     formatting of electrode files between bids-formatted
%                     bair data and non-bids-formatted SOM data)
% specs.plotmesh    = flag to plot meshes with atlases: yes/no
% specs.plotlabel   = flag to plot electrode labels on mesh: yes/no                    
%
% OUTPUT: a struct containing the following fields for two probablistic
% atlases of visual brain regions (wang2015 and benson14)
% out.area_names    = list of all areas in the atlas
% out.area_count    = number of nodes with matches in area_names
% out.elec_labels   = names of electrodes with matches in visual regions
% out.area_labels   = names of areas with matched electrodes
% node_indices      = indices of the matched electrodes; the benson14 atlas
%                     additionally contains eccentricity, polar angle and
%                     sigma estimates for these nodes (probabilistic)
%   
% EXAMPLE INPUT:
% specs.pID      = 'ny648';
% specs.elecFile = '/Volumes/server/Projects/BAIR/Data/BIDS/visual/sub-ny648/ses-NYUECOG01/ieeg/sub-ny648_ses-NYUECOG01_electrodes.tsv';
% specs.fsDir    = '/Volumes/server/Freesurfer_subjects/som648';
% specs.thresh   = []; % default is 20 mm
% specs.patientPool = 'BAIR';
% specs.plotmesh  = 'yes'; % plot meshes with atlases for each subject: yes or no
% specs.plotlabel = 'yes'; % plot electrode labels on mesh: yes or no
%
% EXAMPLE OUTPOUT:
% out.benson14_varea
%   struct with fields:
%       area_names: {'V1'  'V2'  'V3'  'hV4'  'VO1'  'VO2'  'LO1'  'LO2'  'TO1'  'TO2'  'V3b'  'V3a'}
%       area_count: [0 2 0 0 0 1 0 0 2 3 1 4]
%      elec_labels: {'G05'  'G06'  'G07'  'G08'  'G15'  'G16'  'G24'  'IO01'  'MO01'  'MO02'  'MO03'  'MO04'  'DPIL06'}
%      area_labels: {'V3a'  'V3b'  'TO1'  'TO1'  'TO2'  'TO2'  'TO2'  'VO2'  'V2'  'V2'  'V3a'  'V3a'  'V3a'}
%     node_indices: [157670 157047 155681 155786 163093 164068 169918 186636 143611 143710 144737 147429 154577]
%       node_eccen: [1.09 8.29 1.29 0.29 12 0.12 0.08 0.24 4.85 3.85 2.96 10.23 12]
%       node_angle: [19 142.41 143.17 159.98 69.79 60.39 113.53 99.7 140.68 92.95 162.34 113.24 176.61]
%       node_sigma: [1.39 3.96 1.94 0.57 6 0.21 0.14 0.49 0.97 0.8 2.06 4.66 6]

pID    = specs.pID;
thresh = specs.thresh;

if isempty(thresh)
    thresh = inf;
end

plotlabel = specs.plotlabel;
plotmesh  = specs.plotmesh;

% RUN
disp(['running patient ' num2str(pID)]);

% Read electrode coordinate file from BAIR/ECOG directory
elec_file = specs.elecFile;

switch specs.patientPool
    case 'BAIR'
        if exist(elec_file, 'file')

            % Prefer to use readtable for tsv files because it doesn't require
            % knowing the order of the columns beforehand (as textscan does). 
            E = readtable(elec_file, 'FileType', 'text');
            elec_labels = E.name;
            if iscell(E.x)
                for ii = 1:length(elec_labels)
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
        D = dir([patientDir num2str(pID) '/*coor_T1*.txt']);
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
surf_file_rh = [specs.fsDir '/surf/rh.pial'];
surf_file_lh = [specs.fsDir '/surf/lh.pial'];
% FOR TESTING COLORMAPS (if we want to plot inflated/sphere versions, need
% to adjust center coordinates which default to zero in freesurfer
%surf_file_rh = [specs.fsDir '/surf/rh.sphere'];
%surf_file_lh = [specs.fsDir '/surf/lh.sphere'];

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
atlasNames = {'wang2015_atlas', 'benson14_varea', 'benson14_eccen', 'benson14_angle', 'benson14_sigma', 'template_areas'};

% Output
out.patientID = pID;
    
for a = 1:length(atlasNames)
    
    currentAtlas = atlasNames{a};
    
    % Get atlases for this subject
    atlas_file_rh = [specs.fsDir '/surf/rh.' currentAtlas '.mgz'];
    atlas_file_lh = [specs.fsDir '/surf/lh.' currentAtlas '.mgz'];
    if exist(atlas_file_rh, 'file') && exist(atlas_file_lh, 'file')
        [atlas_rh] = load_mgh(atlas_file_rh);
        [atlas_lh] = load_mgh(atlas_file_lh);
    else
        disp(['no annotations found for ' currentAtlas]);
        out.(currentAtlas) = [];
        continue
    end
    
    % Match nearest nodes to atlas labels
    atlas = [squeeze(atlas_rh);squeeze(atlas_lh)]; % concatenate hemis   
    atlas_elec = atlas(indices);
    [elec_indices] = find(atlas_elec);
    elec_labels_found = elec_labels(elec_indices);
    node_indices = indices(elec_indices);
    
    % Clip Benson atlases (results in better colormap scaling)
    switch currentAtlas
        case 'benson14_eccen'
            atlas(atlas>12) = 12.00;
        case 'benson14_sigma'
            atlas(atlas>6) = 6.00;
        case 'benson14_angle'
            %area_cmap_rh = A.cmap(:,6:8);
    end
    
    if ~isempty(elec_labels_found)

        % Get names / colormaps associated with each atlas
        switch currentAtlas

            case 'wang2015_atlas' % Wang atlas

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
                    case 'yes'
                        
                        A = load(['colormap_' currentAtlas]);
                        disp(['loading colormap_' currentAtlas]);
                        area_cmap = A.cmap(:,2:4);

                        %atlas_range = range(atlas);
                        %cmap_index = round(linspace(1,length(area_cmap),atlas_range));
                        %area_cmap = area_cmap(cmap_index,:);  
                    
                    otherwise
                        % do nothing
                end
        end
        
        % Get area labels and specs for matched nodes
        switch currentAtlas
            case {'wang2015_atlas', 'benson14_varea', 'template_areas'}

                area_count = zeros(1,length(area_labels));
                for i = 1:length(elec_labels_found)
                    area_count(atlas_elec(elec_indices(i))) = area_count(atlas_elec(elec_indices(i))) + 1;
                end
                area_labels_found = area_labels(round(atlas_elec(elec_indices)));

                out.(currentAtlas).area_count   = area_count;
                out.(currentAtlas).elec_labels  = elec_labels_found';
                out.(currentAtlas).area_labels  = area_labels_found;
                out.(currentAtlas).node_indices = node_indices;

            case 'benson14_eccen'
                out.benson14_varea.node_eccen = round(atlas(out.benson14_varea.node_indices),2)';

            case 'benson14_angle' 
                out.benson14_varea.node_angle = round(atlas(out.benson14_varea.node_indices),2)';

            case 'benson14_sigma'
                out.benson14_varea.node_sigma = round(atlas(out.benson14_varea.node_indices),2)';
        end
    else
        disp(['no electrodes in ' currentAtlas]);
        out.(currentAtlas) = [];
    end
        
    % Plot
    switch plotmesh
        case 'yes'
            figure('Name', [num2str(pID) ' ' currentAtlas]); hold on;
            
            plot_electrodes(elec_xyz, [1 1 1]*0.2,2);
            plot_electrodes(elec_xyz(elec_indices,:), [0 0 0],2);

            t_r = trimesh(faces_r+1, vertices_r(:,1), vertices_r(:,2), vertices_r(:,3), atlas_rh, 'FaceColor', 'flat'); 
            t_r.LineStyle = 'none';
            axis equal; hold on;
            t_l = trimesh(faces_l+1, vertices_l(:,1), vertices_l(:,2), vertices_l(:,3), atlas_lh, 'FaceColor', 'flat'); 
            t_l.LineStyle = 'none';             
            cmap = [[1 1 1]*.7; area_cmap];
            colormap(cmap); 
            caxis([0 max(atlas)]); %length(area_cmap)]);

            plot_electrodes(vertices(indices,:), [1 1 1]*0.8, 1);
            plot_electrodes(vertices(indices(elec_indices),:), [1 1 1], 1);
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
            % Do not plot
    end
end


% Print to window how many maps were found, and
% make a count per area (across hemispheres)

% node indices to atlases
atlasNames = {'wang2015_atlas', 'benson14_varea', 'template_areas'};

for a = 1:length(atlasNames)
    currentAtlas = atlasNames{a};
    
    if ~isempty(out.(currentAtlas)) && ~isempty(out.(currentAtlas).elec_labels)
        
        disp(['found ' num2str(length(out.(currentAtlas).elec_labels)) ' electrodes in ' currentAtlas ': '])
            
        for i = 1:length(out.(currentAtlas).elec_labels)            
            switch currentAtlas
                case {'wang2015_atlas', 'template_areas'};
                    disp([out.(currentAtlas).elec_labels{i} ' in area ' out.(currentAtlas).area_labels{i}]);
                case 'benson14_varea'
                    disp([out.(currentAtlas).elec_labels{i} ' in area ' out.(currentAtlas).area_labels{i} ...
                        ' with eccen = ' num2str(out.benson14_varea.node_eccen(i)) ', angle = ' num2str(out.benson14_varea.node_angle(i)) ', sigma = ' num2str(out.benson14_varea.node_sigma(i))]);
            end
        end
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

end