
function [out] = electrode_to_nearest_node(specs, varargin)

% This function matches a list of electrode locations in an ECoG patient to
% the nearest node in their T1s freesurfer pial surface reconstruction,
% determines which of those nodes fall within a set of visual regions
% specified by several probabilistic atlases, and extracts node
% information. Note that electrode coordinates should be specified in the
% T1 volume space that was used to obtain the freesurfer reconstruction.
%
% REQUIRED INPUT: a struct containing the following fields:
% specs.pID         = patient ID 
%
% OPTIONAL additional fields: 
% specs.atlasNames  = list of atlases to match electrode locations with.
%                     Current options are (default = all):
%                                         {'wang2015_atlas', ...
%                                          'benson14_varea', ...
%                                          'benson14_eccen', ...
%                                          'benson14_angle', ...
%                                          'benson14_sigma', ...
%                                          'template_areas'};
%                     NOTE: including benson14_varea is required to
%                     be able to obtain benson14_eccen, angle and sigma
% specs.plotmesh    = flag to plot meshes with atlases: 'left', 'right',
%                     'both', 'none';
% specs.plotelecs   = flag to plot electrodes on mesh: 'yes', 'no'
% specs.plotlabel   = flag to plot electrode labels on mesh: 'yes', 'no'
% specs.plotcbar    = flag to plot colorbar for the atlas: 'yes', 'no'
% specs.thresh      = maximum allowed distance between electrode and node,
%                     in mm (if empty, thresh is infinite, meaning that the
%                     electrode can be infinitely far) - default empty
% specs.fsDir       = freesurfer directory of patient (include full path);
%                     should contain the wang and benson atlases that can
%                     be obtained through the nben/neuropythy docker >
%                     default '/Volumes/server/Freesurfer_subjects/som
% specs.patientPool = should be either 'BAIR' or 'SOM' (may be removed
%                     later but necessary now to deal with differences in
%                     formatting of electrode files between bids-formatted
%                     and non-bids-formatted data) - default SOM
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
% specs.pID           = '668'; % patient ID number
% specs.plotmesh      = 'left'; % left, right, both, or none
% specs.plotelecs     = 'no'; % yes or no
% out = electrode_to_nearest_node(specs);
% 
% EXAMPLE OUTPOUT:
% out.benson14_varea
%   struct with fields:
%       area_names: {'V1'  'V2'  'V3'  'hV4'  'VO1'  'VO2'  'LO1'  'LO2'  'TO1'  'TO2'  'V3b'  'V3a'}
%       area_count: [0 0 0 0 0 0 0 4 2 0 1 0]
%      elec_labels: {'LTO01'  'LTO02'  'LTO03'  'RTO01'  'RTO02'  'DLLT08'  'DLLT09'}
%      area_labels: {'LO2'  'LO2'  'TO1'  'V3b'  'LO2'  'TO1'  'LO2'}
%     node_indices: [154623 159762 168662 11390 20685 177358 169529]
%       node_eccen: [0.39 0.21 0.11 2.42 9.45 5.82 27.62]
%       node_angle: [106.68 178.32 177.33 136.13 80.88 0 64.06]
%       node_sigma: [0.75 0.61 0.33 1.87 7.52 8.16 21.11]
      
if ~isfield(specs, 'atlasNames') || isempty(specs.atlasNames)
    specs.atlasNames = {'wang2015_atlas', 'benson14_varea', 'benson14_eccen', 'benson14_angle', 'benson14_sigma', 'template_areas'};
end

if ~isfield(specs, 'patientPool') || isempty(specs.patientPool)
    specs.patientPool = 'SOM'; 
end

if ~isfield(specs, 'fsDir') || isempty(specs.fsDir)
    switch specs.patientPool
        case 'SOM'
            specs.fsDir = fullfile(filesep, 'Volumes', 'server', 'Freesurfer_subjects', sprintf('som%s',specs.pID)); 
        case 'BAIR'
            specs.fsDir = fullfile('/Volumes/server/Freesurfer_subjects/', num2str(specs.pID)); 
    end
end

if ~isfield(specs, 'thresh') || isempty(specs.thresh)
    specs.thresh = inf;
end

if ~isfield(specs, 'plotmesh') || isempty(specs.plotmesh)
    specs.plotmesh = 'both';
end

if ~isfield(specs, 'plotelecs') || isempty(specs.plotelecs)
    switch specs.plotmesh
        case 'none'
            specs.plotelecs = 'no';
        otherwise
            specs.plotelecs = 'yes';
    end
end

if ~isfield(specs, 'plotlabel') || isempty(specs.plotlabel)
    specs.plotlabel = 'yes';
end

if ~isfield(specs, 'plotcbar') || isempty(specs.plotcbar)
    specs.plotcbar = 'yes';
end

plotmesh  = specs.plotmesh;
plotelecs = specs.plotelecs;
plotlabel = specs.plotlabel;
plotcbar  = specs.plotcbar;

 % Do we have a patient ID?
if ~isfield(specs, 'pID') || isempty(specs.pID)
    disp('please specify a patient ID');
    return
else
    disp(['running patient ' num2str(specs.pID)]);
end
        
% Read electrode coordinate file from BAIR or RAW directory
switch specs.patientPool
    case 'BAIR'
        
        % For BIDS formatted data, patient ID is specified in dataDir
        patientDir = varargin{1};
        
        D = dir([patientDir '/*electrodes*.tsv']);
        if ~isempty(D)
           
            elec_file = D(1).name;
            disp(['reading ' elec_file]); 

            % Prefer to use readtable for tsv files because it doesn't require
            % knowing the order of the columns beforehand (as textscan does). 
            
            E = readtable(fullfile(D(1).folder, elec_file), 'FileType', 'text');
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
        
        % Check whether this patient is in the RAW BAIR directory, if not find it in SoM
        patientDir = '/Volumes/server/Projects/BAIR/Data/Raw/ECoG/';
        if ~isdir([patientDir num2str(specs.pID)])
            patientDir = '/Volumes/server/Projects/BAIR/Data/Raw/ECoG/SoM/';
            if ~isdir([patientDir num2str(specs.pID)])
                disp('patient directory not found - exiting');
                out = [];
                return
            end
        end
        
        % Read electrode coordinate file from Raw/SoM directory
        D = dir([patientDir num2str(specs.pID) '/*coor_T1*.txt']);
        if ~isempty(D)
            elec_file = D(1).name;
            disp(['reading ' elec_file]); % if there are multiple coor_T1 files for separate hemisphere, D(1) will always be the full list
            fid = fopen([patientDir num2str(specs.pID) '/' elec_file]); E = textscan(fid, '%s%f%f%f%s'); fclose(fid);
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
keep_idx = find(bestSqDist<specs.thresh);
indices = indices(keep_idx);
elec_xyz = elec_xyz(keep_idx,:);

% Output
out.patientID = specs.pID;


for a = 1:length(specs.atlasNames)
    
    currentAtlas = specs.atlasNames{a};
    
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
                    case {'both', 'left', 'right'}
                    
                        A = load(['colormap_' currentAtlas]);
                        disp(['loading colormap_' currentAtlas]);
                        
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
                atlasUnits = 'visual area';

            case 'benson14_eccen'
                out.benson14_varea.node_eccen = round(atlas(out.benson14_varea.node_indices),2)';
                atlasUnits = 'degrees';
                % bin eccentricity values?
%                atlas_rh = ceil(atlas_rh);
%                atlas_lh = ceil(atlas_lh);
%                atlas_rh(atlas_rh > 20) = 20;

            case 'benson14_angle' 
                out.benson14_varea.node_angle = round(atlas(out.benson14_varea.node_indices),2)';
                atlasUnits = 'degrees';
                % atlas(atlas>6) = 6.00;

            case 'benson14_sigma'
                out.benson14_varea.node_sigma = round(atlas(out.benson14_varea.node_indices),2)';
                atlasUnits = 'degrees';
        end
    else
        disp(['no electrodes in ' currentAtlas]);
        out.(currentAtlas) = [];
    end
        
    % Plot
    switch plotmesh
        
        case 'none'
            % Do nothing
            
        otherwise
            
            % Plot mesh
            fig = figure('Name', [num2str(specs.pID) ' ' currentAtlas]); hold on;
            
            switch plotmesh
                case 'both'          
                    
                    t_r = trimesh(faces_r+1, vertices_r(:,1), vertices_r(:,2), vertices_r(:,3), atlas_rh, 'FaceColor', 'flat'); 
                    t_r.LineStyle = 'none';
                    hold on;
                    t_l = trimesh(faces_l+1, vertices_l(:,1), vertices_l(:,2), vertices_l(:,3), atlas_lh, 'FaceColor', 'flat'); 
                    t_l.LineStyle = 'none';             
                    
                case 'left'

                    t_l = trimesh(faces_l+1, vertices_l(:,1), vertices_l(:,2), vertices_l(:,3), atlas_lh, 'FaceColor', 'flat'); 
                    t_l.LineStyle = 'none';             
                
                case 'right'
                    
                    % Plot mesh
                    t_r = trimesh(faces_r+1, vertices_r(:,1), vertices_r(:,2), vertices_r(:,3), atlas_rh, 'FaceColor', 'flat'); 
                    t_r.LineStyle = 'none';
                    
            end
            
            axis equal; hold on;
            cmap = [[1 1 1]*.7; area_cmap];
            colormap(gcf,cmap);
                        
             % Clip Benson atlases (results in better colormap scaling)
            switch currentAtlas
                case 'benson14_eccen'
                    caxis([0 20]);
                case 'benson14_sigma'
                    caxis([0 10]);
                otherwise
                    caxis([0 max(atlas)]); %length(area_cmap)]);
            end
            
            % Add colorbar?
            switch plotcbar
                case 'yes'
                    cb = colorbar;
                    cb.FontSize = 18;
                    cb.Color = [0 0 0];
                    cb.Label.String = atlasUnits;
                    cb.Position = [0.92 0.25 0.03 0.5];
                    switch currentAtlas
                        case {'wang2015_atlas', 'benson14_varea', 'template_areas'}
                            cb.Ticks = 0:1:length(area_labels);
                            cb.TickLabels = ['none', area_labels];
                    end   
            end
    end
    
    switch plotelecs
        
        case 'yes'            
            if ~exist('fig','var')
                fig = figure('Name', [num2str(specs.pID) ' ' currentAtlas]); hold on;
            end

            switch plotmesh
                case 'left'
                    electoplot = find(elec_xyz(:,1) < 0);
                    elec_indices = intersect(elec_indices,electoplot);
                    elec_plotindex = electoplot;
                case 'right'
                    electoplot = find(elec_xyz(:,1) >= 0);
                    elec_indices = intersect(elec_indices,electoplot);
                    elec_plotindex = electoplot;
                otherwise
                    elec_plotindex = 1:size(elec_xyz,1);
            end
            
            % Plot electrodes
            plot_electrodes(elec_xyz(elec_plotindex,:), [1 1 1]*0.2,2);
            plot_electrodes(elec_xyz(elec_indices,:), [0 0 0],2);

            % Plot matched nodes
            plot_electrodes(vertices(indices(elec_plotindex),:), [1 1 1]*0.8, 1);
            plot_electrodes(vertices(indices(elec_indices),:), [1 1 1], 1);
            
            switch plotlabel
                case 'yes'
                    for i = 1:size(elec_xyz(elec_plotindex,:),1)
                        [x, y, z] = adjust_elec_label(elec_xyz(elec_plotindex(i),:),2);
                        text('Position',[x y z],'String',elec_labels(elec_plotindex(i),:),'Color','w','VerticalAlignment','top');
                    end
            end
            
        otherwise
            % do nothing
    end
           
    if exist('fig','var')
        
        % Set view parameters
        axis off; set(gcf, 'color','white','InvertHardCopy', 'off');
        view(0,0);
        material dull;

        h=light; lightangle(h,  45, 45); lighting gouraud;
        h=light; lightangle(h, -45, 45); lighting gouraud;
        h=light; lightangle(h, -45, -90); lighting gouraud;

        set(gcf, 'Position', [150 100 1500 1250]);
        axis tight
    end
end

% Print to window how many maps were found, and
% make a count per area (across hemispheres)

% Check which of these atlases we ran
atlasNames = intersect({'wang2015_atlas', 'benson14_varea', 'template_areas'}, specs.atlasNames);

% In case not all benson maps were run, fill fields to prevent error below
if max(contains(atlasNames, 'benson14_varea')) &&  ~isempty(out.benson14_varea)
    if ~isfield(out.benson14_varea, 'node_eccen'); out.benson14_varea.node_eccen = nan(1,length(out.benson14_varea.node_indices));end
    if ~isfield(out.benson14_varea, 'node_angle'); out.benson14_varea.node_angle = nan(1,length(out.benson14_varea.node_indices));end
    if ~isfield(out.benson14_varea, 'node_sigma'); out.benson14_varea.node_sigma = nan(1,length(out.benson14_varea.node_indices));end
end      
    
for a = 1:length(atlasNames)
    currentAtlas = atlasNames{a};
    
    if ~isempty(out.(currentAtlas)) && ~isempty(out.(currentAtlas).elec_labels)
        
        disp(['found ' num2str(length(out.(currentAtlas).elec_labels)) ' electrodes in ' currentAtlas ': '])
            
        for i = 1:length(out.(currentAtlas).elec_labels)            
            switch currentAtlas
                case {'wang2015_atlas', 'template_areas'}
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


%% PATH TO NOAH COLORMAPS:
% '/Volumes/server/Projects/HCP/analysis/images'


% alignment to fs average
% find nearest node ID in pial surface of indiv.subj
% find coordinate of this node ID in indiv.subj lh.sphere.reg
% find nearest node ID with closest coordinate in fsaverage lh.sphere
% that node ID gives you coordinates on lh.pial of fsaverage

end