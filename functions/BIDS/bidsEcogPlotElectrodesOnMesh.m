function varargout = bidsEcogPlotElectrodesOnMesh(projectDir, subject, session, atlasName, varargin)
% bidsEcogPlotElectrodesOnMesh(projectDir, subject, session, atlasName, [thresh], [surfaceType], [specs])
% 
% Plots iEEG electrode positions from a BIDS directory onto on a 3D brain
% mesh of the pial surface reconstruction in derivatives/freesurfer, along
% with an retinotopic atlas (if specified and present in derivatives).
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
%                           - glasser16_atlas
%                           - none (default: plot mesh without atlas)
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
%     specs:
%     specs.plotmesh    = flag to plot meshes with atlases: 'auto'(default),
%                         'left', 'right', 'both' (if 'auto', plot hemisphere
%                         where electrodes locate)
%     specs.plotelecs   = flag to plot electrodes on mesh: 'yes', 'no'
%     specs.plotlabel   = flag to plot electrode labels on mesh: 'yes', 'no'
%     specs.plotmatchednodes = flag to plotmatched nodes in mesh: 'yes', 'no'
%     specs.labelsize   = font size of electrode labels (default = [])
%     specs.plotcbar    = flag to plot colorbar for the atlas: 'yes', 'no'
%     specs.adjustLRdistance = flag to adjust the distance between left and right surface: 'auto', 'yes', 'no'
%     specs.face_alpha  = transparency value for mesh (default 1);
%     specs.plotelecrad = radius of electrodes. If empty, will use the
%                         electrode size in the electrode table as radius for
%                         each electrode. If electrode size is not provided,
%                         will automatically scale all electrodes to 2, or to
%                         user specified scalar (specified here).
%     specs.plotnoderad = radius of matched node (default = 1)
%     specs.view        = camera angles of plotting (default = [0,0])
%     specs.areaName    = cell-array of visual area names, which are
%                         plotted when full probability map is specified
%                         (default = 'all')
% 
% Example 1
% This example plots electrode positions for subject p10 on the normalized 
% wang full probability atlas and benson atlas
%     projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual_ecog_recoded'; 
%     subject           = 'p10';
%     atlasName         = {'wang15_fplbl_norm_5', 'benson14_varea'};
%     specs = [];
%     specs.plotelecrad = 1;
% h = bidsEcogPlotElectrodesOnMesh(projectDir, subject, [], atlasName, [], specs);
% 
% Example 2
% This example plots electrode positions for subject p10 on the wang full
% probability atlas for 'V1', 'V2', and 'V3' with marging detailed visual areas
%     projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual_ecog_recoded'; 
%     subject           = 'p10';
%     atlasName         = {'smry_wang15_fplbl'};
%     specs = [];
%     specs.plotelecrad = 1;
%     specs.areaName    = {'V1','V2','V3'};
% h = bidsEcogPlotElectrodesOnMesh(projectDir, subject, [], atlasName, [], specs);
% 
% Example 3
% This example plots electrode positions for subject p10 on the wang atlas on the white matter
%     projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual_ecog_recoded'; 
%     subject           = 'p10';
%     atlasName         = 'wang15_mplbl';
%     surfaceType       = 'white';
%     specs = [];
%     specs.plotelecs   = 'yes';
%     specs.plotelecrad = 1;
% h = bidsEcogPlotElectrodesOnMesh(projectDir, subject, [], atlasName, [], surfaceType, specs);
% 
% See also, bidsSpecifySessions, bidsEcogReadElectrodeFile,
%           bidsEcogReadSurfFile, bidsEcogReadAtlasFile, bidsEcogGetMatchedAtlas,
%           bidsEcogMatchElectrodesToAtlas

% KY, BAIR 2022

% <projectDir>
if ~exist('projectDir', 'var') || isempty(projectDir), error('projectDir not defined'); end    

% <subject>
if ~exist('subject', 'var') || isempty(subject), error('subject not defined'); end

% <subject>
if ~exist('session', 'var') || isempty(session)
	[sessions] = bidsSpecifySessions(projectDir, subject);
    idx = find(contains(lower(sessions), {'ecog', 'iemu'}));
    if isempty(idx), error('no ECOG sessions found for this subject');end
    session = sessions{idx(1)};
end

% <atlas>
if ~exist('atlasName', 'var') || isempty(atlasName), atlasName = 'nane'; end
if ~iscell(atlasName), atlasName = {atlasName}; end

% <optionals>
thresh      = [];
surfaceType = [];
specs       = [];
if numel(varargin)>2
    thresh      = varargin{1};
    surfaceType = varargin{2};
    specs       = varargin{3};
elseif numel(varargin)>1
    if ischar(varargin{1})||isstring(varargin{1})||isstring(varargin{1})
        surfaceType = varargin{1};
        specs       = varargin{2};
    else
        thresh      = varargin{1};
        if isstruct(varargin{2})
            specs       = varargin{2};
        else
            surfaceType = varargin{2};
        end
    end
elseif numel(varargin)>0
    if ischar(varargin{1})||isstring(varargin{1})
        surfaceType = varargin{1};
    elseif isstruct(varargin{1})
        specs       = varargin{1};
    else
        thresh      = varargin{1};
    end
end

% <surfaceOption>
surfaceBase        = [];
surfaceTransSmooth = [];
if isstruct(surfaceType)
    if isfield(surfaceType,'surfaceTransSmooth')
        surfaceTransSmooth = surfaceType.surfaceTransSmooth;
    end
    if isfield(surfaceType,'surfaceBase')
        surfaceBase        = surfaceType.surfaceBase;
    end
    if isfield(surfaceType,'surfaceType')
        surfaceType        = surfaceType.surfaceType;
    else
        surfaceType        = [];
    end
end

% <thresh>
if isempty(thresh), thresh = inf; end

% <surfaceType>
if isempty(surfaceType)
    if isfield(specs, 'surf') && ~isempty(specs.surf)
        surfaceType = specs.surf;
    elseif isfield(specs, 'surfaceType') && ~isempty(specs.surfaceType)
        surfaceType = specs.surfaceType;
    else
        % surfaceType should be 'pial' to match electrode locations in BAIR project
        surfaceType = 'pial';
    end
end
if isempty(surfaceBase) && isfield(specs, 'surfaceBase')
    surfaceBase        = specs.surfaceBase;
end
if isempty(surfaceTransSmooth) && isfield(specs, 'surfaceTransSmooth')
    surfaceTransSmooth = specs.surfaceTransSmooth;
end

%% Specs
if ~isfield(specs, 'plotmesh') || isempty(specs.plotmesh)
    specs.plotmesh = 'auto';
end

if ~isfield(specs, 'plotelecs') || isempty(specs.plotelecs)
    switch specs.plotmesh
        case 'none'
            specs.plotelecs = 'no';
            specs.plotcbar = 'no';
        otherwise
            switch surfaceType
                case 'pial'
                    specs.plotelecs = 'yes';
                otherwise
                    specs.plotelecs = 'no';
            end
    end
end
if ~isfield(specs, 'plotelecrad')
    specs.plotelecrad = [];
end

if ~isfield(specs, 'plotlabel') || isempty(specs.plotlabel)
    specs.plotlabel = 'no';
end

if ~isfield(specs, 'labelsize')
    specs.labelsize = [];
end

if ~isfield(specs, 'plotmatchednodes') || isempty(specs.plotmatchednodes)
    specs.plotmatchednodes = 'no';
end

if ~isfield(specs, 'plotnoderad')
    specs.plotnoderad = [];
end

if ~isfield(specs, 'plotcbar') || isempty(specs.plotcbar)
    specs.plotcbar = 'yes';
end

if ~isfield(specs, 'face_alpha') || isempty(specs.face_alpha)
    specs.face_alpha = 1;
end

if ~isfield(specs, 'view') || isempty(specs.view)
    specs.view = [0,0];
end

if ~isfield(specs, 'areaName') || isempty(specs.areaName)
    specs.areaName = 'all';
end

if ~isfield(specs, 'adjustLRdistance') || isempty(specs.adjustLRdistance)
    specs.adjustLRdistance = 'auto';
end

plotmesh         = specs.plotmesh;
face_alpha       = specs.face_alpha;
plotelecs        = strcmpi(specs.plotelecs,'yes');
plotlabel        = strcmpi(specs.plotlabel,'yes');
plotcbar         = strcmpi(specs.plotcbar,'yes');
plotmatchednodes = strcmpi(specs.plotmatchednodes,'yes');
adjustLRdistance = specs.adjustLRdistance;
plotelecrad      = specs.plotelecrad;
plotnoderad      = specs.plotnoderad;
viewanlge        = specs.view;
if isempty(specs.labelsize)
    labelsize = {};
else
    labelsize = {'FontSize',specs.labelsize};
end
areaName        = specs.areaName;
if ~iscell(areaName),  areaName = {areaName};  end

%% Read in files and matched ECoG electrodes and MRI atlas

[~,~,~,~,atlasName] = interpretAtlasNames(atlasName);
[matched_atlas_vals, electrode_table, matched_vertices, keep_idx, ~, ...
    elec_xyz, atlases_r, atlases_l, vertices_r, faces_r, vertices_l, faces_l, atlasName] = ...
    bidsEcogGetMatchedAtlas(projectDir, subject, session, atlasName, thresh, surfaceType, surfaceBase, surfaceTransSmooth);
elec_labels = electrode_table.name(keep_idx);
elec_size = electrode_table.size(keep_idx);

%% Prepare figures

% Estimate electrode size
if isempty(plotelecrad)
    if exist('elec_size', 'var')
        elec_radii = elec_size;
    else
        elec_radii = ones(size(elec_xyz,1),1)*2;
    end
else
    elec_radii = ones(size(elec_xyz,1),1)*plotelecrad;
end
if isempty(plotnoderad)
    node_radii = elec_radii*0.75;
else
    node_radii = ones(size(elec_xyz,1),1)*plotnoderad;
end

% Estimate appropriate hemisphere
elec_isr = ismember(electrode_table.hemisphere,'R');
if strcmpi(plotmesh,'auto')
    plotmesh = 'both';
    if plotelecs
        elecbias = sum(elec_isr) ./ size(matched_vertices,1);
        if     elecbias(1) == 1,  plotmesh = 'right';
        elseif elecbias(1) == 0,  plotmesh = 'left';
        end
    end
end

% Adjust surface location
surfdist = min(vertices_r(:,1)) - max(vertices_l(:,1));
surfwdth = diff(minmax(vertcat(vertices_l(:,1),vertices_r(:,1))'));
if strcmpi(adjustLRdistance,'auto')
    adjustLRdistance = surfdist < -surfwdth/10 | surfdist > surfwdth/10;
else
    adjustLRdistance = strcmpi(adjustLRdistance,'yes');
end
if adjustLRdistance
    isodist  = [round((-surfdist + surfwdth/20)./2) 0 0];
    elec_xyz   = elec_xyz + (elec_isr-~elec_isr).*isodist;
    matched_vertices ...
               = matched_vertices + (elec_isr-~elec_isr).*isodist;
    vertices_r = vertices_r + isodist;
    vertices_l = vertices_l - isodist;
end

%% Make Figures

h = gobjects(0);
figidx = 1;
for a = 1:length(atlasName)
    currentAtlas = atlasName{a};
    fprintf('[%s] Plotting atlas %s\n',mfilename,currentAtlas);
    
    % Look up area names
    [area_labels, area_cmap, units] = getAtlasLabels(currentAtlas);
    
    % Get atlas values
      matched_atlas = matched_atlas_vals{a};
      atlas_rh      = atlases_r{a};
      atlas_lh      = atlases_l{a};
      
    % check full probability map
    isfplbl = ~iscolumn(matched_atlas);
     
    % Convert Noah's porlar anlge to analyzePRF porlar angle
    switch currentAtlas
        case {'benson14_angle'}
            atlas_rh = bensonng2PRFang(atlas_rh,'right');
            atlas_lh = bensonng2PRFang(atlas_lh,'left');
    end
    
    % Match electrodes with atlas
    if isfplbl              % atlas is vertices x M (FPM)
        ncells = size(matched_atlas);  ncells(1) = 1;
        elec_indices      = cell(ncells);
        for i = 1:prod(ncells(2:end))
            elec_indices{i}      = find(matched_atlas(:,i));
        end
    else                    % atlas is vertices x 1 (MPM or others)
        [elec_indices] = find(matched_atlas);
    end
    
    % Get names / colormaps associated with each atlas
    if isempty(area_labels)
        [area_cmap, units, col_range] = elecrend_ReadColormap(currentAtlas);
    end
  
    % Plot
    for i = 1:size(matched_atlas,2)
        if ~isfplbl || any(ismember([area_labels(i), {'all'}], areaName))
            if isfplbl
                    figName = sprintf('%s %s %s',num2str(subject),currentAtlas,area_labels{i});
            elseif  ismember(currentAtlas, {'none'})
                    figName = sprintf('%s',num2str(subject));
            else
                    figName = sprintf('%s %s',num2str(subject),currentAtlas);
            end

            switch plotmesh

                case 'none'
                    % Do nothing

                otherwise

                    % Plot mesh
                    h(figidx) = figure('Name', figName); hold on;

                    switch plotmesh
                        case 'both'                                                 
                            plot_mesh(faces_r, vertices_r, atlas_rh(:,i), area_cmap, face_alpha);
                            plot_mesh(faces_l, vertices_l, atlas_lh(:,i), area_cmap, face_alpha);
                        case 'left'
                            plot_mesh(faces_l, vertices_l, atlas_lh(:,i), area_cmap, face_alpha);
                        case 'right'
                            plot_mesh(faces_r, vertices_r, atlas_rh(:,i), area_cmap, face_alpha);   
                    end

                     % Clip Benson atlases (results in better colormap scaling)
                    switch units
                        case {'area','parcellation'}
                            caxis([0 size(area_cmap,1)+1]); 
                        case 'probability'
                            caxis([0 1]); 
                        case 'degrees'
                            caxis(col_range);
                        otherwise
                            caxis([0 max(matched_atlas(:,i))]); %length(area_cmap)]);
                    end

                    % Add colorbar
                    if plotcbar && ~isempty(area_cmap)
                        cb = colorbar;
                        cb.FontSize = 18;
                        cb.Color = [0 0 0];
                        cb.Label.String = units;
                        cb.Position = [0.92 0.25 0.03 0.5];
                        switch units
                            case {'area','parcellation'}
                                dTick = diff(caxis)/(size(area_cmap,1)+1);
                                cb.Ticks = (min(caxis)+dTick/2):dTick:(max(caxis));
                                cb.TickLabels = ['none', area_labels];
                        end   
                    end
            end

            if plotelecs

                if ~ishandle(h(figidx))
                    h(figidx) = figure('Name', figName); hold on;
                end

                % Get hemi electrodes
                elec_plotindex = getElecinHemi(elec_xyz, plotmesh, 10);
                if iscell(elec_indices)
                    elec_selindices = intersect(elec_indices{i},elec_plotindex);
                else
                    elec_selindices = intersect(elec_indices,elec_plotindex);
                end

                % Plot electrodes
                plot_electrodes(elec_xyz(elec_plotindex,:), [1 1 1]*0.2,elec_radii(elec_plotindex));
                plot_electrodes(elec_xyz(elec_selindices,:), [0 0 0],elec_radii(elec_selindices));

                % Plot matched nodes
                if plotmatchednodes
                    plot_electrodes(matched_vertices(elec_plotindex,:), [1 1 1]*0.8, node_radii(elec_plotindex));
                    plot_electrodes(matched_vertices(elec_selindices,:), [1 1 1], node_radii(elec_plotindex));
                end

                if plotlabel
                    for j = 1:size(elec_xyz(elec_plotindex,:),1)
                        [x, y, z] = adjust_elec_label(elec_xyz(elec_plotindex(j),:),elec_radii(elec_plotindex(j)).*1.3,plotmesh);
                        text('Position',[x y z],'String',elec_labels(elec_plotindex(j),:),'Color','w','VerticalAlignment','top',labelsize{:});
                    end
                end
            end

            if numel(h) == figidx
                % Set view parameters
                set_view(gcf,viewanlge);
                figidx = figidx + 1;
            end
        end
    end
end
if nargout > 0
    varargout{1} = h;
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% SUBROUTINES%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_electrodes(xyz, color, radius)   
    [x1, y1, z1] = sphere;
    %if ~exist('radius', 'var'), radius = 2; end
    if ~exist('radius', 'var') || isempty(radius), radius = ones(size(xyz,1),1)*2; end 
    if isscalar(radius), radius = ones(size(xyz,1),1)*radius;end
    %x1 = radius * x1;
    %y1 = radius * y1;
    %z1 = radius * z1;
    for jj = 1:size(xyz,1)
        s = surface(xyz(jj,1)+x1*radius(jj), ...
                    xyz(jj,2)+y1*radius(jj), ...
                    xyz(jj,3)+z1*radius(jj), ...
                    'FaceColor', color, 'EdgeColor', 'none', 'AmbientStrength',0.7);
        s.CDataMapping = 'direct';
    end
end

function plot_mesh(faces, vertices, atlas, area_cmap, alpha)
    t = trimesh(faces+1, vertices(:,1), vertices(:,2), vertices(:,3), atlas, 'FaceColor', 'flat', 'FaceAlpha', alpha); 
    t.LineStyle = 'none';    
	axis equal; hold on;
    %if cmapHasZeroIndex
        cmap = [[1 1 1]*.7; area_cmap];
    %else
    %    cmap = area_cmap;
    %end
    colormap(gcf,cmap);    
end

function set_view(gcf,viewangle)
    if ~exist('viewangle','var') || isempty(viewangle)
        viewangle = [0,0];
    end
    axis off; set(gcf, 'color','white','InvertHardCopy', 'off');
    view(viewangle);
    material dull;
    h=light; lightangle(h,  45, 45); lighting gouraud;
    h=light; lightangle(h, -45, 45); lighting gouraud;
    h=light; lightangle(h, -45, -90); lighting gouraud;
    set(gcf, 'Position', [150 100 1500 1250]);
    axis tight    
end

function [x, y, z] = adjust_elec_label(xyz,radius,plotmesh)
    if ~exist('radius','var') || isempty(radius)
        radius = 3;
    end
    if ~exist('plotmesh','var') || isempty(plotmesh)
        plotmesh = 'unknown';
    end
    
    switch plotmesh
        case 'right'
            xbias = 5;
        case 'left'
            xbias = -5;
        otherwise
            xbias = 0;
    end
        
    if xyz(1)>xbias
        x = xyz(1)+radius;
    else
        x = xyz(1)-radius;
    end
    if xyz(3)>0
        z = xyz(3)+radius;
    else
        z = xyz(3)-radius;
    end
    if abs(xyz(1))<abs(xbias)
        y = xyz(2)-radius;
    else
        y = xyz(2);
    end
end

function [ang] = bensonng2PRFang(ang, hemi)
    % convert benson angle to pRF angle
    
    ang(ang == 0) = nan;              % zero indicates no data
    switch lower(hemi)
        case {'left','l'}
            ang = mod(90 - ang - 1, 360) + 1; % Right Visual Field (avoiding 0)
        case {'right','r'}
            ang = mod(90 + ang - 1, 360) + 1; % Left Visual Field (avoiding 0)
    end
    ang(isnan(ang)) = 0;              % put back zero    
    
end

function [elec_plotindex] = getElecinHemi(elec_xyz, hemi, allowance)
    % get electrode index in hemisphere
    
    if nargin < 3 || isempty(allowance), allowance = 0; end
    
    switch lower(hemi)
        case {'left','l'}
            elec_plotindex = find(elec_xyz(:,1) <= allowance);
        case {'right','r'}
            elec_plotindex = find(elec_xyz(:,1) >= -allowance);
        otherwise
            elec_plotindex = find(elec_xyz(:,1) <= allowance |...
                                  elec_xyz(:,1) >= -allowance);
    end
    
end