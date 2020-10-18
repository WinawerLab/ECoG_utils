function [electrode_table] = bidsEcogMatchElectrodesToAtlas(projectDir, subject, session, atlasName, thresh, printSummary)
% Matches electrodes to visual atlas, using electrode coordinates specified
% in the bids metadata file in electrodes.tsv and freesurfer surface nodes
% in derivates/freesurfer. Area names will be added to the electrode_table.
%
% [electrode_table] = bidsEcogMatchElectrodesToAtlas(projectDir, ...
%   subject, session, atlas, [thresh], [printSummary]);
%
% Notes:
% Atlases are expected to be located in freesurfer's surf directory and
% formatted as .mgz files; they can be obtained using Noah Benson's
% retinotopy atlas docker, see run_docker_retinotopyatlas.sh in
% https://github.com/BAIRR01/BAIRanalysis/tree/master/scripts
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
%                           - glasser 
%     thresh:           maximum allowed distance between electrode and node,
%                       in mm (if empty, thresh is infinite, meaning that the
%                       electrode can be infinitely far) 
%                           default empty
%     printSummary:     prints number of matched electrodes and matched
%                       areas (if applicable) to command windown
%                           default true 
%
%
% Example 1
% This example matches electrode positions for subject som726 to the wang
% maximum probability atlas
%     projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual'; 
%     subject           = 'som726';
%     atlasName         = 'wang15_mplbl';
% electrode_table = bidsEcogMatchElectrodesToAtlas(projectDir, subject, [], atlasName);
% 
% Example 2
% This example matches electrode positions for subject som726 to the wang
% full probability atlas and benson atlas
%     projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual'; 
%     subject           = 'som726';
%     atlasName         = {'wang15_fplbl', 'benson14_varea'};
% electrode_table = bidsEcogMatchElectrodesToAtlas(projectDir, subject, [], atlasName);
%
% See also bidsSpecifySessions.m bidsEcogReadElectrodeFile.m
%          bidsEcogReadSurfFile.m bidsEcogReadAtlasFile.m
% Makes use of 'nearpoints' function from vistasoft

% IG, BAIR 2020

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
if ~exist('atlasName', 'var') || isempty(atlasName), error('atlasName not defined. See help for list of options.'); end

% <thresh>
if ~exist('thresh', 'var') || isempty(thresh), thresh = inf; end

% <printSummary>
if ~exist('printSummary', 'var') || isempty(printSummary), printSummary = true; end

%% Read in files

% Read in electrode cooordinates
[electrode_table] = bidsEcogReadElectrodeFile(projectDir, subject, session);

% Read surface reconstructions from BIDS derivatives 
[vertices_r, ~, vertices_l, ~] = bidsEcogReadSurfFile(projectDir, subject);

% Read in atlases
if ~iscell(atlasName), atlasName = {atlasName}; end
for a = 1:length(atlasName)
    % Atlases can have different dimensions, so we're reading them in as
    % cell arrays
    [atlases_r{a}, atlases_l{a}] = bidsEcogReadAtlasFile(projectDir, subject, atlasName{a});
end

%% Match electrodes to nodes

% Prepare data for matching
elec_labels = electrode_table.name;
elec_xyz    = [electrode_table.x electrode_table.y electrode_table.z];
vertices    = [vertices_r;vertices_l];

% TEMP: fixing an issue with sub-chaam's coordinates, should not be in here but
% should be fixed in the files themselves 
if contains(subject,{'umcuchaam', 'chaam'})
    % subtract effect of cropping by freesurfer
    elec_xyz(:,1) = elec_xyz(:,1)-3.4490;
    elec_xyz(:,2) = elec_xyz(:,2)-34.6040;
    elec_xyz(:,3) = elec_xyz(:,3)+6.2660;
end

% Match the electrode xyz with the nodes in the surfaces; find nearest
[indices, bestSqDist] = nearpoints(elec_xyz', vertices'); % function from vistasoft

% Ignore electrodes that are more than [thresh] away from any surface node
keep_idx = find(bestSqDist < thresh);
indices = indices(keep_idx);
elec_xyz = elec_xyz(keep_idx,:);

%% Match nodes to atlases 

% Match nearest nodes to atlas labels
matched_atlas_vals = cell(length(atlasName),1);
for a = 1:length(atlasName)
    
    atlas_rh = squeeze(atlases_r{a});   
    atlas_lh = squeeze(atlases_l{a});
    if ~iscolumn(atlas_rh)
        atlas_rh = permute(atlas_rh,[ndims(atlas_rh), 1:(ndims(atlas_rh)-1)]); % align with rows
        atlas_lh = permute(atlas_lh,[ndims(atlas_lh), 1:(ndims(atlas_lh)-1)]); % align with rows
    end
    atlas = cat(1,atlas_rh,atlas_lh); % concatenate hemis   
    matched_atlas_vals{a} = atlas(indices,:);
end

% Add matched nodes to electrode table
electrode_table.node_x = vertices(indices,1);
electrode_table.node_y = vertices(indices,2);
electrode_table.node_z = vertices(indices,3);
electrode_table.node_idx = indices';

% Add atlas info to electrode table
for a = 1:length(atlasName)
    
    vals = matched_atlas_vals{a};
    
    % Add column to electrode table with the atlas values for the matched
    % nodes, or the area names corresponding to those atlas values
    
     % Look up area names
    [area_labels] = getAtlasLabels(atlasName{a});
    
    if isempty(area_labels)
        % This atlas has no area labels; just add the node values themselves
        electrode_table.(atlasName{a}) = vals;
    else        
        if iscolumn(vals) % if single column, replace node values with area name
            matched_areas = cell(size(vals));
            idx = vals>0;
            matched_areas(idx) = area_labels(vals(idx));
            matched_areas(~idx) = {'none'};
            electrode_table.(atlasName{a}) = matched_areas;
        else % otherwise, add separate column for each area with atlas vals
            ncols = size(vals,2);
            for ii = 1:ncols
                colName = sprintf('%s_%s', atlasName{a}, area_labels{ii});
                electrode_table.(colName) = vals(:,ii);
            end
        end
    end
end

if printSummary
    
    idx = contains(electrode_table.Properties.VariableNames, atlasName);
    colNames = electrode_table.Properties.VariableNames(idx);
	
    % Print a summary for each column
    for ii = 1:length(colNames)
        vals = electrode_table.(colNames{ii});
        if isnumeric(vals)
            nMatched = length(find(vals>0));
            fprintf('[%s] %s: %d matched electrodes out of %d electrodes total \n', ...
                mfilename, colNames{ii}, nMatched, height(electrode_table));
        else
            elec_idx = find(~contains(vals, 'none'));
            nMatched = length(elec_idx);
            fprintf('[%s] %s: %d matched electrodes out of %d electrodes total: \n', ...
                mfilename, colNames{ii}, nMatched, height(electrode_table));
            % Print a match for each electrode for with atlas area name
            for jj = 1:length(elec_idx)
                fprintf('[%s] %s in area %s\n', mfilename, ...
                electrode_table.name{elec_idx(jj)}, vals{elec_idx(jj)})
            end
           
        end
    end
end
