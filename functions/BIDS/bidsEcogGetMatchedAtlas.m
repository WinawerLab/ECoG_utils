function  [matched_atlas_vals, electrode_table, matched_vertices, keep_idx, indices, ...
           elec_xyz, atlases_r, atlases_l, vertices_r, faces_r, vertices_l, faces_l, atlasName] = ...
             bidsEcogGetMatchedAtlas(projectDir, subject, session, atlasName, thresh, surfaceType, surfaceBase, surfaceTransSmooth)
% [matched_atlas_vals, electrode_table, matched_vertices, keep_idx, indices, ...
%  elec_xyz, atlases_r, atlases_l, vertices_r, faces_r, vertices_l, faces_l, atlasName] = ...
%   bidsEcogGetMatchedAtlas(projectDir, subject, session, atlasName, thresh, surfaceType, [surfaceBase, surfaceTransSmooth])
% 
% Read electrodes and atlas and match them in space.
% 
% Input
%     projectDir:       path where the BIDS project lies (string)
%     subject:          BIDS subject name (string, all lower case)
%     session:          BIDS session name (string, all lower case)
%                          if empty: first session with 'ecog' in the name
%     atlasName:        name of atlas to match electrode locations with. Can
%                       also be a list of names formatted as a cell array.
%                       current options are:
%                           - wang15_mplbl
%                           - wang15_fplbl
%                           - benson14_varea
%                           - benson14_angle
%                           - benson14_eccen
%                           - benson14_sigma
%                           - benson20_mplbl
%                           - benson20_fplbl
%                           - glasser 
%                           - none (plot mesh without atlas)
%                       additional options are: (see also interpretAtlasNames)
%                           - wang15_fplbl_norm_*
%                           - benson20_fplbl_norm_*
%                               apply normalization with *% threshold
%                           - smry_wang15_...
%                           - smry_benson14_...
%                               merge some visual areas for simplification
%     thresh:           maximum allowed distance between electrode and node,
%                       in mm (if empty, thresh is infinite, meaning that the
%                       electrode can be infinitely far) 
%                           default empty
%     surfaceType:      which MRI surface to read 
%                           default pial
% 
% Supplementary Input
%     surfaceBase:          default pial
%     surfaceTransSmooth:   default 10 [mm]
% 
% 
% See also, bidsEcogMatchElectrodesToAtlas, bidsEcogPlotElectrodesOnMesh
% 
% Makes use of 'nearpoints' function from vistasoft
         
% K.Yuasa, IG, BAIR 2022;
        

% <thresh>
if ~exist('thresh', 'var') || isempty(thresh), thresh = inf; end
if ~exist('transrefdist', 'var'), surfaceTransSmooth = []; end

% <surfaceType>
if ~exist('surfaceType','var') || isempty(surfaceType)
    surfaceType = 'pial';   % surface for output
end
if ~exist('surfaceBase','var') || isempty(surfaceBase)
    surfaceBase = 'pial';   % surface to corregister with electrodes
end

%% Read in files

% Read in electrode cooordinates
[electrode_table, elec_xyz] = bidsEcogReadElectrodeFile(projectDir, subject, session);

% Read surface reconstructions from BIDS derivatives 
[vertices_r, faces_r, vertices_l, faces_l] = bidsEcogReadSurfFile(projectDir, subject, surfaceType);
if ~strcmpi(surfaceType,surfaceBase)
    [basevertices_r, ~, basevertices_l] = bidsEcogReadSurfFile(projectDir, subject, surfaceBase);
else
    basevertices_r = vertices_r;    basevertices_l = vertices_l;
end

% Interpret atlas name
if ~iscell(atlasName), atlasName = {atlasName}; end
[atlasFileName, isnorm, normthresh,issmry,atlasName] = interpretAtlasNames(atlasName);

% Read in atlases
atlases_r = cell(1,length(atlasFileName)); atlases_l = cell(1,length(atlasFileName));
for a = 1:length(atlasFileName)
    % Atlases can have different dimensions, so we're reading them in as
    % cell arrays
    if strcmpi(atlasFileName{a},'none')
        atlases_r{a} = zeros(1,1,size(vertices_r,1));
        atlases_l{a} = zeros(1,1,size(vertices_l,1));
    elseif find(ismember(atlasFileName,atlasFileName{a}),1) < a
        loadeda = find(ismember(atlasFileName,atlasFileName{a}),1);
        atlases_r(a) = atlases_r(loadeda);
        atlases_l(a) = atlases_l(loadeda);
    else
        [atlases_r{a}, atlases_l{a}] = bidsEcogReadAtlasFile(projectDir, subject, atlasFileName{a});
    end
end

%% Reshape atlas values

% Match nearest nodes to atlas labels
for a = 1:length(atlasName)
    
    atlas_rh = squeeze(atlases_r{a});   
    atlas_lh = squeeze(atlases_l{a});
    % Apply normalization
    if isnorm(a) % full probability map normalized
        atlas_rh = normalizefplbl(atlas_rh, normthresh(a));
        atlas_lh = normalizefplbl(atlas_lh, normthresh(a));
    end
    if ~iscolumn(atlas_rh)
        atlas_rh = permute(atlas_rh,[ndims(atlas_rh), 1:(ndims(atlas_rh)-1)]); % align with rows
        atlas_lh = permute(atlas_lh,[ndims(atlas_lh), 1:(ndims(atlas_lh)-1)]); % align with rows
    end
    % Apply area marging
    [~, ~, ~, area_match] = getAtlasLabels(atlasFileName{a},isnorm(a),issmry(a));
    if ~iscolumn(atlas_rh)      % i.e. full probability map
        maxarea   = max(area_match(:));
        area_idx  = zeros(1,maxarea);
        for ilabel=reshape(unique(area_match),1,[])
            mergearea = find(area_match==ilabel);
            area_idx(ilabel) = mergearea(1);
            atlas_rh(:,area_idx(ilabel)) = sum(atlas_rh(:,mergearea),2);
            atlas_lh(:,area_idx(ilabel)) = sum(atlas_lh(:,mergearea),2);
        end
        atlas_rh = atlas_rh(:,area_idx);
        atlas_lh = atlas_lh(:,area_idx);
    else
        maxarea   = length(area_match);
        for ilabel=1:maxarea
            atlas_rh(atlas_rh==ilabel) = area_match(ilabel) + maxarea;
            atlas_lh(atlas_lh==ilabel) = area_match(ilabel) + maxarea;
        end
        atlas_rh(atlas_rh>maxarea) = atlas_rh(atlas_rh>maxarea) - maxarea;
        atlas_lh(atlas_lh>maxarea) = atlas_lh(atlas_lh>maxarea) - maxarea;
    end
    % Put back
    atlases_r{a}          = atlas_rh;
    atlases_l{a}          = atlas_lh;
    
end

%% Match electrodes to nodes

% Isolate left and right hemispheres not to be associated to the opposite
if ismember('hemisphere',fieldnames(summary(electrode_table))) && ...
        all(ismember(electrode_table.hemisphere,{'L','R'}))
    elec_isr = ismember(electrode_table.hemisphere,'R');
    surfdist = min(basevertices_r(:,1)) - max(basevertices_l(:,1));
    surfwdth = diff(minmax(vertcat(basevertices_l(:,1),basevertices_r(:,1))'));
    isodist  = [round(max(-surfdist,0) + surfwdth/20) 0 0];
else
    elec_isr = false(height(electrode_table),1);
    isodist  = [0 0 0];
end
elec_xyz_iso   = elec_xyz + elec_isr.*isodist;
basevertices   = [basevertices_r + isodist; basevertices_l];
            
% Match the electrode xyz with the nodes in the surfaces; find nearest
[indices, bestSqDist] = nearpoints(elec_xyz_iso', basevertices'); % function from vistasoft

% Transform electrode location
if ~strcmpi(surfaceType,surfaceBase) && ismember('inflated',{surfaceType,surfaceBase})
    if  ~ismember(surfaceBase,{'smoothwm','inflated'})   % surfaceType = 'inflated'
        [basevertices_r, ~, basevertices_l] = bidsEcogReadSurfFile(projectDir, subject, 'smoothwm');
        basevertices   = [basevertices_r + isodist; basevertices_l];
        [tmpindices] = nearpoints(elec_xyz_iso', basevertices'); % function from vistasoft
        tmpvertices_r = vertices_r;     tmpvertices_l = vertices_l;
    elseif ~ismember(surfaceType,{'smoothwm','inflated'})  % surfaceBase = 'inflated'
        [tmpvertices_r, ~, tmpvertices_l] = bidsEcogReadSurfFile(projectDir, subject, 'smoothwm');
        tmpindices = indices;
    end
    vertices     = [tmpvertices_r;tmpvertices_l];
    basevertices = [basevertices_r; basevertices_l];
        
    [sulc_r, sulc_l] = bidsEcogReadCurvFile(projectDir, subject, 'sulc');
    sulc = [sulc_r;sulc_l];
    
    elec_xyz = transform_elec_xyz(elec_xyz,tmpindices,basevertices,vertices,size(basevertices_r,1),sulc,surfaceTransSmooth);
    
end
vertices   = [vertices_r;vertices_l];

% Add LR information
if ~ismember('hemisphere',fieldnames(summary(electrode_table)))
    elec_isr = indices <= size(vertices_r,1);
    electrode_table.hemisphere(~elec_isr) = {'L'};
    electrode_table.hemisphere(elec_isr)  = {'R'};
end

% Ignore electrodes that are more than [thresh] away from any surface node
keep_idx = find(bestSqDist < thresh);
indices  = indices(keep_idx);
elec_xyz = elec_xyz(keep_idx,:);

%% Match nodes to atlases 

% Match nearest nodes to atlas labels
matched_vertices   = vertices(indices,:);
matched_atlas_vals = cell(length(atlasName),1);
for a = 1:length(atlasName)
    
    atlas_rh = atlases_r{a};   
    atlas_lh = atlases_l{a};
    atlas   = cat(1,atlas_rh,atlas_lh); % concatenate hemis   
    matched_atlas_vals{a} = atlas(indices,:);
    
end
% 
% % Add matched nodes to electrode table
% electrode_table.node_x(keep_idx) = matched_vertices(:,1);
% electrode_table.node_y(keep_idx) = matched_vertices(:,2);
% electrode_table.node_z(keep_idx) = matched_vertices(:,3);
% electrode_table.node_idx(keep_idx) = indices';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% SUBROUTINES%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [atlas] = normalizefplbl(atlas, normthresh)
    % assign each node with probability > normthresh to an area
    tmp = zeros(1,size(atlas,2));
    ind = sum(atlas,1) > max(normthresh/100,0);
    [~,tmp(:,ind)] = max(atlas(:,ind));
    atlas = tmp;
end

function [elec_xyz] = transform_elec_xyz(elec_xyz,indices,basevertices,vertices,Nrhemi,sulc,transsmooth)
    
    if ~exist('transsmooth','var')||isempty(transsmooth), transsmooth = 10; end
    if ~exist('sulc','var'), sulc = []; end
    transsmooth = max(transsmooth,0);
    
    % make cluster to compute transform matrix
    isrhemi    = false(size(vertices,1),1); isrhemi(1:Nrhemi) = true;
    elec_isr = indices <= Nrhemi;
    for el = 1:numel(indices)
        vertdist    = [vecnorm(basevertices - basevertices(indices(el),:),2,2),...
                       vecnorm(vertices - vertices(indices(el),:),2,2)];
        indices_sel = all(vertdist <= transsmooth,2) & isrhemi==elec_isr(el);
        if ~isempty(sulc)   % segregate sulcus
            indices_sel = indices_sel & (sulc<-1)==(sulc(indices(el))<-1);
        end
        transmat    = preptrns(basevertices(indices_sel,:))\preptrns(vertices(indices_sel,:));
        elec_xyznew = preptrns(elec_xyz(el,:))*transmat;
        
        elec_xyz(el,:) = elec_xyznew(:,1:size(elec_xyz,2));
    end

end

function [xyzmat] = preptrns(xyzmat)
    xyzmat = cat(2,xyzmat,ones(size(xyzmat,1),1));
end