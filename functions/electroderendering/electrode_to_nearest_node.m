function [out] = electrode_to_nearest_node(specs, BIDSformatted)

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
% specs.rootDir     = directory that has the patient electrode positions, 
%                     default: /Volumes/server/Projects/BAIR/Data/BIDS. If
%                     not found in default, looks in Raw ECoG data folder
%                     (/Volumes/server/Projects/BAIR/Data/Raw/ECoG).
% specs.fsDir       = freesurfer directory of patient (include full path);
%                     default: /Volumes/server/Freesurfer_subjects/som
%                     for plotting atlases, the freesurfer directory should 
%                     contain the wang and benson atlases that can be 
%                     obtained through the nben/neuropythy docker
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
% BIDSformatted     = flag to indicate if the data is in BIDS
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
% specs.pID           = 'som668'; % patient ID 
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Check inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(specs, 'rootDir') || isempty(specs.rootDir)
    specs.rootDir = fullfile(filesep, 'Volumes', 'server', 'Projects', 'BAIR', 'Data');
    if ~exist('BIDSformatted', 'var')
        if contains(specs.pID, {'som', 'umcu', 'ny'})
            % this is a BIDS formatted dataset; look in the BIDS folder
            BIDSformatted = 1;
        else
            % this is not a BIDS formatted dataset; look in the Raw folder 
            BIDSformatted = 0;
        end
    end
end

if ~isfield(specs, 'fsDir') || isempty(specs.fsDir)
%     specs.fsDir = fullfile(filesep, 'Volumes', 'server', 'Freesurfer_subjects'); 
    specs.fsDir = fullfile(filesep, 'Volumes', 'server', 'Projects','BAIR','Data','BIDS','visual','derivatives','freesurfer');
end

if ~isfield(specs, 'atlasNames') || isempty(specs.atlasNames)
    specs.atlasNames = {'wang2015_atlas','wang15_mplbl', 'benson14_varea', 'benson14_eccen', 'benson14_angle', 'benson14_sigma', 'template_areas'};
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
            specs.plotcbar = 'no';
        otherwise
            specs.plotelecs = 'yes';
    end
end
if ~isfield(specs, 'plotelecrad') || isempty(specs.plotelecrad)
    specs.plotelecrad = 2;
end

if ~isfield(specs, 'plotlabel') || isempty(specs.plotlabel)
    specs.plotlabel = 'yes';
end

if ~isfield(specs, 'plotmatchednodes') || isempty(specs.plotmatchednodes)
    specs.plotmatchednodes = 'no';
end
if ~isfield(specs, 'plotnoderad') || isempty(specs.plotnoderad)
    specs.plotnoderad = 1;
end

if ~isfield(specs, 'plotcbar') || isempty(specs.plotcbar)
    specs.plotcbar = 'yes';
end

 % Do we have a patient ID?
if ~isfield(specs, 'pID') || isempty(specs.pID)
    fprintf('[%s] Please specify a patient ID\n',mfilename);
    return
else
    fprintf('[%s] Running electrode_to_node for patient %s...\n',mfilename,num2str(specs.pID));
end

plotmesh         = specs.plotmesh;
plotelecs        = specs.plotelecs;
plotlabel        = specs.plotlabel;
plotcbar         = specs.plotcbar;
plotmatchednodes = specs.plotmatchednodes;
plotelecrad      = specs.plotelecrad;
plotnoderad      = specs.plotnoderad;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Read in files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if BIDSformatted
    
    specs.dataDir = fullfile(specs.rootDir, 'BIDS'); 

    % Check which projectfolder this patient is located
    patientDir = fullfile(specs.dataDir, 'visual', sprintf('sub-%s', specs.pID));
    if ~isfolder(patientDir)
        patientDir = fullfile(specs.dataDir, 'motor', sprintf('sub-%s', specs.pID));
    end
    if ~isfolder(patientDir)
        patientDir = fullfile(specs.dataDir, 'tactile', sprintf('sub-%s', specs.pID));
    end
    patientDir = fullfile(specs.dataDir, 'pilots', 'visualPilot', sprintf('sub-%s', specs.pID));
    if ~isfolder(patientDir)
        patientDir = fullfile(specs.dataDir, 'pilots', 'motorPilot', sprintf('sub-%s', specs.pID));
    end
    if ~isfolder(patientDir)
        patientDir = fullfile(specs.dataDir, 'pilots', 'tactilePilot', sprintf('sub-%s', specs.pID));
    end
    if ~isfolder(patientDir)
        patientDir = fullfile(specs.rootDir, sprintf('sub-%s', specs.pID));
    end
    if ~isfolder(patientDir)
        fprintf('[%s] Could not find patient data folder. Please check paths. \n',mfilename);
        out = [];
        return
    end

    % Find an ECoG session directory
    sessionList = dir(fullfile(patientDir, sprintf('ses*ecog*')));
    if isempty(sessionList)
        sessionList = dir(fullfile(patientDir, sprintf('ses*ECOG*')));
    end
    if isempty(sessionList)
        fprintf('[%s] Could not find ECoG sessions for this patient. Please check paths. \n',mfilename);
        out = [];
        return
    else
        
        % Read electrode coordinate file
        D = dir(fullfile(patientDir, sessionList(1).name, 'ieeg', '*electrodes.tsv'));          
        if isempty(D)
            fprintf('[%s] Electrode coordinate file not found in %s - exiting [%s]. \n', ...
                mfilename,fullfile(patientDir, sessionList(1).name, mfilename));
            out = [];
            return
        else
            elec_file = D(1).name;
            elec_folder = D(1).folder;
            fprintf('[%s] Reading %s...\n',mfilename,fullfile(elec_folder,elec_file));

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
            if contains(specs.pID,{'umcuchaam', 'chaam'})
                % subtract effect of cropping by freesurfer
                elec_xyz(:,1) = elec_xyz(:,1)-3.4490;
                elec_xyz(:,2) = elec_xyz(:,2)-34.6040;
                elec_xyz(:,3) = elec_xyz(:,3)+6.2660;
                % add to1 and to2 probability atlases
                specs.atlasNames  = [specs.atlasNames {'to1_percent', 'to2_percent'}];
                % add glasser atlas
                specs.atlasNames  = [specs.atlasNames 'glasser16_atlas'];
            end
        end
    end

else
    
    specs.dataDir = fullfile(rootDir, 'Raw', 'ECoG'); 
    patientDir = fullfile(specs.dataDir, specs.pID);

    % Check whether this patient is in the RAW BAIR directory, if not find it in SoM
    if ~isfolder(patientDir)
        specs.dataDir = '/Volumes/server/Projects/BAIR/Data/Raw/ECoG/SoM/';
        patientDir = fullfile(specs.dataDir, specs.pID);
        if ~isfolder(patientDir)
            fprintf('[%s] Patient data directory not found - exiting [%s]. \n',mfilename,mfilename);
            out = [];
            return
        end
    end

    % Read electrode coordinate file from Raw/SoM directory
    D = dir(fullfile(patientDir, '/*coor_T1*.txt'));
    if ~isempty(D)
        elec_file = D(1).name;% if there are multiple coor_T1 files for separate hemisphere, D(1) will always be the full list
        elec_folder = D(1).folder;
        fprintf('[%s] Reading %s...\n',mfilename,fullfile(elec_folder,elec_file));
        fid = fopen(fullfile(patientDir,elec_file)); E = textscan(fid, '%s%f%f%f%s'); fclose(fid);
        elec_xyz = [E{2} E{3} E{4}]; 
        elec_labels = E{1};
        elec_types = unique(E{5});
        fprintf('[%s] Types of electrodes found: %s\n',mfilename,[elec_types{:}]);
    else
        fprintf('[%s] Electrode coordinate file not found - exiting [%s] (NYU: is the server mapped?)\n',mfilename,mfilename);
        out = [];
        return
    end
    
    % Add 'som' to name in order to be able to read freesurfer directory
    specs.pID = sprintf('som%s',specs.pID);
end

% Read surface reconstructions from Freesurfer_subjects directory
if contains(specs.fsDir, 'BIDS') || BIDSformatted; specs.pID = sprintf('sub-%s', specs.pID);end
surf_file_rh = fullfile(specs.fsDir, specs.pID, 'surf', 'rh.pial');
surf_file_lh = fullfile(specs.fsDir, specs.pID, 'surf', 'lh.pial');

if exist(surf_file_rh, 'file') && exist(surf_file_lh, 'file')
    fprintf('[%s] Reading Freesurfer pial surface file %s ...\n',mfilename, surf_file_rh);
    [vertices_r, faces_r] = read_surf(surf_file_rh);
    fprintf('[%s] Reading Freesurfer pial surface file %s ...\n',mfilename, surf_file_lh);
    [vertices_l, faces_l] = read_surf(surf_file_lh);
    vertices = [vertices_r;vertices_l];
else
    fprintf('[%s] Freesurfer surfaces not found - exiting [%s] \n',mfilename,mfilename);
    out = [];
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Match elecs to nodes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Match the electrode xyz with the nodes in the surfaces; find nearest
[indices, bestSqDist] = nearpoints(elec_xyz', vertices'); % function from vistasoft

% Ignore electrodes that are more than [thresh] away from any surface node
keep_idx = find(bestSqDist<specs.thresh);
indices = indices(keep_idx);
elec_xyz = elec_xyz(keep_idx,:);

% Initalize output
out.patientID = specs.pID;

%global cmapHasZeroIndex 
%cmapHasZeroIndex = true;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Generate figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FIRST, plot electrodes on brain without any atlas %%%%%%%%%%%%%%%%%%%%

switch plotmesh
    case 'none'
        % do nothing
    otherwise
        fig = figure('Name', specs.pID); hold on;
end

% Plot pial surface as mesh
switch plotmesh
    case 'both'                                                 
        plot_mesh(faces_l, vertices_l, ones([1 1 size(vertices_l,1)]), []);
        plot_mesh(faces_r, vertices_r, ones([1 1 size(vertices_r,1)]), []);
    case 'left'
        plot_mesh(faces_l, vertices_l, ones([1 1 size(vertices_l,1)]), []);
    case 'right'
        plot_mesh(faces_r, vertices_r, ones([1 1 size(vertices_r,1)]), []);   
end

% Add the electrodes
switch plotelecs        
    case 'yes'            
        if ~exist('fig','var')
            fig = figure('Name', [num2str(specs.pID) ' ' currentAtlas]); hold on;
        end
        
        [elec_indices] = find(indices);

        switch plotmesh
            case 'left'
                electoplot = find(elec_xyz(:,1) <= 10);
                elec_indices = intersect(elec_indices,electoplot);
                elec_plotindex = electoplot;
            case 'right'
                electoplot = find(elec_xyz(:,1) >= -10);
                elec_indices = intersect(elec_indices,electoplot);
                elec_plotindex = electoplot;
            otherwise
                elec_plotindex = 1:size(elec_xyz,1);
        end

        % Plot electrodes
        plot_electrodes(elec_xyz(elec_plotindex,:), [1 1 1]*0.2,2);
        plot_electrodes(elec_xyz(elec_indices,:), [0 0 0],2);
        
        % Plot matched nodes
        switch plotmatchednodes
            case 'yes'
                plot_electrodes(vertices(indices(elec_plotindex),:), [1 1 1]*0.8, 1);
                plot_electrodes(vertices(indices(elec_indices),:), [1 1 1], 1);
        end
            
        % Add electrode labels
        switch plotlabel
            case 'yes'
                for i = 1:size(elec_xyz(elec_plotindex,:),1)
                    [x, y, z] = adjust_elec_label(elec_xyz(elec_plotindex(i),:),3);
                    text('Position',[x y z],'String',elec_labels(elec_plotindex(i),:),'Color','w','VerticalAlignment','top');
                end
        end
end

% Set view parameters
if exist('fig','var')
    set_view(gcf)
end

% SECOND, match electrodes with each atlas %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make figure if requested %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for a = 1:length(specs.atlasNames)
    
    currentAtlas = specs.atlasNames{a};  
    % Get atlases for this subject
    %fullfile(specs.fsDir, specs.pID, 'surf', 'rh.pial');
    atlas_file_rh = fullfile(specs.fsDir, specs.pID, 'surf', sprintf('rh.%s.mgz', currentAtlas));
    atlas_file_lh = fullfile(specs.fsDir, specs.pID, 'surf', sprintf('lh.%s.mgz', currentAtlas));
    if exist(atlas_file_rh, 'file') && exist(atlas_file_lh, 'file')
        [atlas_rh] = load_mgh(atlas_file_rh);
        [atlas_lh] = load_mgh(atlas_file_lh);
    else
        fprintf('[%s] No annotations found for %s \n',mfilename, currentAtlas);
        out.(currentAtlas) = [];
        continue % move on to the next atlas
    end
    
    % Match nearest nodes to atlas labels
    atlas_rh = squeeze(atlas_rh);   atlas_lh = squeeze(atlas_lh);
    if ~iscolumn(atlas_rh)  % 1st dimension should be vertices 
        atlas_rh = permute(atlas_rh,[ndims(atlas_rh), 1:(ndims(atlas_rh)-1)]); % align with rows
        atlas_lh = permute(atlas_lh,[ndims(atlas_lh), 1:(ndims(atlas_lh)-1)]); % align with rows
    end
    atlas = cat(1,atlas_rh,atlas_lh); % concatenate hemis   
    atlas_elec = atlas(indices,:);
    if iscolumn(atlas_elec) % atlas_elec is vertices x 1 (MPM)
        [elec_indices] = find(atlas_elec);
        elec_labels_found = elec_labels(elec_indices);
        node_indices = indices(elec_indices);
    else                    % atlas_elec is vertices x M (FPM)
        ncells = size(atlas_elec);  ncells(1) = 1;
        elec_indices      = cell(ncells);
        elec_labels_found = cell(ncells);
        node_indices      = cell(ncells);
        for i = 1:prod(ncells(2:end))
            elec_indices{i}      = find(atlas_elec(:,i));
            elec_labels_found{i} = elec_labels(elec_indices{i});
            node_indices{i}      = indices(elec_indices{i});
        end
    end
    
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
            
        case {'wang15_fplbl'} % Wang full probability map

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
            area_cmap   = autumn(64);
            out.(currentAtlas).area_names   = area_labels;

        case {'benson14_varea'} % Noah template: areas 

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

        case {'template_areas','benson20_mplbl'} % old Noah templates & Noah's HCP template

            area_labels = {'V1', 'V2', 'V3'}; 
            area_cmap   = [255 0 0; 255 128 0; 255 255 0]./255;
            out.(currentAtlas).area_names   = area_labels;
            
        case {'benson20_fplbl'} % Noah's HCP full probability map

            area_labels = {'V1', 'V2', 'V3'}; 
            area_cmap   = autumn(64);
            out.(currentAtlas).area_names   = area_labels;
            
        case 'to1_percent' % Wang probability

            area_labels = {'TO1'}; 
            area_cmap = autumn(64);
            out.(currentAtlas).area_names   = area_labels;

        case 'to2_percent' % Wang probability

            area_labels = {'TO2'}; 
            area_cmap = autumn(64);
            out.(currentAtlas).area_names   = area_labels;
            
        case 'glasser16_atlas' % Wang probability

            area_labels = num2cell(1:180);
            %area_cmap = jet(180);
            area_cmap = repmat([1 1 1]*.7, [180 1]);
            area_cmap(2,:) = [1 0 0]; % MST
            area_cmap(23,:) = [1 0.5 0]; % MT
            area_cmap(137,:) = [0 1 0]; % PHT
            area_cmap(138,:) = [0 1 1]; % PH
            area_cmap(139,:) = [1 0.6 0.6]; % TPOJ1
            area_cmap(140,:) = [1 0.8 0.6]; % TPOJ2
            area_cmap(141,:) = [1 1 0.6]; % TPOJ3
            area_cmap(143,:) = [0 0.25 0.2]; %PGp
            area_cmap(156,:) = [0 0.5 1]; %v4t
            area_cmap(157,:) = [0.6 0.6 0]; %FST
            out.(currentAtlas).area_names   = area_labels;
            %cmapHasZeroIndex = true;
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
        
    % Get area labels and specs for matched nodes
    switch currentAtlas
        case {'wang2015_atlas', 'wang15_mplbl', 'benson14_varea','benson20_mplbl', 'template_areas'}

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

        case 'benson14_angle' 
            out.benson14_varea.node_angle = round(atlas(out.benson14_varea.node_indices),2)';
            atlasUnits = 'degrees';
            % atlas(atlas>6) = 6.00;

        case 'benson14_sigma'
            out.benson14_varea.node_sigma = round(atlas(out.benson14_varea.node_indices),2)';
            atlasUnits = 'degrees';
            
        case {'to1_percent', 'to2_percent'}
            out.(currentAtlas).area_count   = length(elec_labels_found);
            out.(currentAtlas).elec_labels  = elec_labels_found';
            out.(currentAtlas).node_indices = node_indices;
            out.(currentAtlas).node_values  = round(atlas(node_indices),2)';
            atlasUnits = 'overlap (%)';
            
        case {'glasser16_atlas'}
            atlasUnits = 'glasser parcellation';
            out.(currentAtlas).area_count   = length(elec_labels_found);
            out.(currentAtlas).elec_labels  = elec_labels_found';
            out.(currentAtlas).node_indices = node_indices;
            out.(currentAtlas).node_values  = round(atlas(node_indices),2)';
            
        case {'wang15_fplbl', 'benson20_fplbl'}
            for i = 1:length(area_labels)
                out.(currentAtlas).area_count(i)   = length(elec_labels_found{i});
                out.(currentAtlas).elec_labels{i}  = elec_labels_found{i}';
                out.(currentAtlas).node_indices{i} = node_indices{i};
                out.(currentAtlas).node_values{i}  = round(atlas(node_indices{i},i),4)';
            end
            atlasUnits = 'overlap';

    end
  
    % Plot
    for i = 1:size(atlas_elec,2)
        switch currentAtlas
            case {'wang15_fplbl', 'benson20_fplbl'}
                posfix = [' ' area_labels{i}];
            otherwise
                posfix = [];
        end
        
        switch plotmesh

            case 'none'
                % Do nothing

            otherwise

                % Plot mesh
                fig = figure('Name', [num2str(specs.pID) ' ' currentAtlas posfix]); hold on;

                switch plotmesh
                    case 'both'                                                 
                        plot_mesh(faces_r, vertices_r, atlas_rh(:,i), area_cmap);
                        plot_mesh(faces_l, vertices_l, atlas_lh(:,i), area_cmap);
                    case 'left'
                        plot_mesh(faces_l, vertices_l, atlas_lh(:,i), area_cmap);
                    case 'right'
                        plot_mesh(faces_r, vertices_r, atlas_rh(:,i), area_cmap);   
                end

                 % Clip Benson atlases (results in better colormap scaling)
                switch currentAtlas
                    case 'benson14_eccen'
                        caxis([0 20]);
                    case 'benson14_sigma'
                        caxis([0 10]);
                    otherwise
                        caxis([0 max(atlas(:,i))]); %length(area_cmap)]);
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
                            case {'wang2015_atlas', 'wang15_mplbl', 'benson14_varea', 'benson20_mplbl', 'template_areas'}
                                dTick = length(area_labels)/(length(area_labels)+1);
                                cb.Ticks = (dTick/2):dTick:length(area_labels);
                                cb.TickLabels = ['none', area_labels];
                        end   
                end
        end

        switch plotelecs

            case 'yes'            
                if ~exist('fig','var')
                    fig = figure('Name', [num2str(specs.pID) ' ' currentAtlas posfix]); hold on;
                end

                switch plotmesh
                    case 'left'
                        electoplot = find(elec_xyz(:,1) <= 10);
                        if iscell(elec_indices)
                          elec_selindices = intersect(elec_indices{i},electoplot);
                        else
                          elec_selindices = intersect(elec_indices,electoplot);
                        end
                        elec_plotindex = electoplot;
                    case 'right'
                        electoplot = find(elec_xyz(:,1) >= -10);
                        if iscell(elec_indices)
                          elec_selindices = intersect(elec_indices{i},electoplot);
                        else
                          elec_selindices = intersect(elec_indices,electoplot);
                        end
                        elec_plotindex = electoplot;
                    otherwise
                        elec_plotindex = 1:size(elec_xyz,1);
                        if iscell(elec_indices)
                          elec_selindices = elec_indices{i};
                        else
                          elec_selindices = elec_indices;
                        end
                end

                % Plot electrodes
                plot_electrodes(elec_xyz(elec_plotindex,:), [1 1 1]*0.2,plotelecrad);
                plot_electrodes(elec_xyz(elec_selindices,:), [0 0 0],plotelecrad);

                % Plot matched nodes
                switch plotmatchednodes
                    case 'yes'
                        plot_electrodes(vertices(indices(elec_plotindex),:), [1 1 1]*0.8, plotnoderad);
                        plot_electrodes(vertices(indices(elec_selindices),:), [1 1 1], plotnoderad);
                end

                switch plotlabel
                    case 'yes'
                        for i = 1:size(elec_xyz(elec_plotindex,:),1)
                            [x, y, z] = adjust_elec_label(elec_xyz(elec_plotindex(i),:),3);
                            text('Position',[x y z],'String',elec_labels(elec_plotindex(i),:),'Color','w','VerticalAlignment','top');
                        end
                end
        end

        if exist('fig','var')
            % Set view parameters
            set_view(gcf)
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Print and Count% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Print to window how many maps were found, and make a count per area (across hemispheres)

% Check which of these atlases we ran
atlasNames = intersect({'wang2015_atlas','wang15_mplbl', 'benson14_varea', 'benson20_mplbl', 'template_areas'}, specs.atlasNames);

% In case not all benson maps were run, fill fields to prevent error below
if max(contains(atlasNames, 'benson14_varea')) &&  ~isempty(out.benson14_varea)
    if ~isfield(out.benson14_varea, 'node_eccen'); out.benson14_varea.node_eccen = nan(1,length(out.benson14_varea.node_indices));end
    if ~isfield(out.benson14_varea, 'node_angle'); out.benson14_varea.node_angle = nan(1,length(out.benson14_varea.node_indices));end
    if ~isfield(out.benson14_varea, 'node_sigma'); out.benson14_varea.node_sigma = nan(1,length(out.benson14_varea.node_indices));end
end      
    
for a = 1:length(atlasNames)
    currentAtlas = atlasNames{a};
    
    if ~isempty(out.(currentAtlas)) && ~isempty(out.(currentAtlas).elec_labels)
        
        fprintf('[%s] Found %d electrodes in %s :\n' ,mfilename, length(out.(currentAtlas).elec_labels), currentAtlas);
            
        for i = 1:length(out.(currentAtlas).elec_labels)            
            switch currentAtlas
                case {'wang2015_atlas', 'wang15_mplbl', 'benson20_mplbl', 'template_areas'}
                    fprintf('[%s] %s in area %s\n', mfilename, out.(currentAtlas).elec_labels{i}, out.(currentAtlas).area_labels{i})
                    %disp([out.(currentAtlas).elec_labels{i} ' in area ' out.(currentAtlas).area_labels{i}]);
                case 'benson14_varea'
                    fprintf('[%s] %s in area %s with eccen = %4.2f, angle = %4.2f, sigma = %4.2f\n', ...
                        mfilename, out.(currentAtlas).elec_labels{i}, out.(currentAtlas).area_labels{i}, ...
                        out.benson14_varea.node_eccen(i), out.benson14_varea.node_angle(i), out.benson14_varea.node_sigma(i));
                    %disp([out.(currentAtlas).elec_labels{i} ' in area ' out.(currentAtlas).area_labels{i} ...
                        %' with eccen = ' num2str(out.benson14_varea.node_eccen(i)) ', angle = ' num2str(out.benson14_varea.node_angle(i)) ', sigma = ' num2str(out.benson14_varea.node_sigma(i))]);
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% SUBROUTINES%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_electrodes(xyz, color, radius)   
    [x, y, z] = sphere;
    if ~exist('radius', 'var'), radius = 3; end
        x = radius * x;
        y = radius * y;
        z = radius * z;
    for ii = 1:size(xyz,1)
       s = surface(xyz(ii,1)+x, xyz(ii,2)+y, xyz(ii,3)+z, 'FaceColor', color, 'EdgeColor', 'none', 'AmbientStrength',0.7);
       s.CDataMapping = 'direct';
    end
end

function plot_mesh(faces, vertices, atlas, area_cmap)
    t = trimesh(faces+1, vertices(:,1), vertices(:,2), vertices(:,3), atlas, 'FaceColor', 'flat'); 
    t.LineStyle = 'none';    
	axis equal; hold on;
    %if cmapHasZeroIndex
        cmap = [[1 1 1]*.7; area_cmap];
    %else
    %    cmap = area_cmap;
    %end
    colormap(gcf,cmap);    
end

function set_view(gcf)
    axis off; set(gcf, 'color','white','InvertHardCopy', 'off');
    view(0,0);
    material dull;
    h=light; lightangle(h,  45, 45); lighting gouraud;
    h=light; lightangle(h, -45, 45); lighting gouraud;
    h=light; lightangle(h, -45, -90); lighting gouraud;
    set(gcf, 'Position', [150 100 1500 1250]);
    axis tight    
end

function [x, y, z] = adjust_elec_label(xyz,radius)
    if ~exist('radius','var')
        radius = 3;
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
