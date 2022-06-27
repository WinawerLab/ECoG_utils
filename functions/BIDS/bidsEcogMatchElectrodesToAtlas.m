function [electrode_table] = bidsEcogMatchElectrodesToAtlas(projectDir, subject, session, atlasName, thresh, surfaceType, printSummary)
% Matches electrodes to visual atlas, using electrode coordinates specified
% in the bids metadata file in electrodes.tsv and freesurfer surface nodes
% in derivates/freesurfer. Area names will be added to the electrode_table.
%
% [electrode_table] = bidsEcogMatchElectrodesToAtlas(projectDir, ...
%   subject, session, atlas, [thresh], [surfaceType], [printSummary]);
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
%                           - benson20_mplbl
%                           - benson20_fplbl
%                           - glasser16_atlas
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
%     printSummary:     prints number of matched electrodes and matched
%                       areas (if applicable) to command windown
%                           default true 
%
%
% Example 1
% This example matches electrode positions for subject p10 to the wang
% maximum probability atlas
%     projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual'; 
%     subject           = 'p10';
%     atlasName         = 'wang15_mplbl';
% electrode_table = bidsEcogMatchElectrodesToAtlas(projectDir, subject, [], atlasName);
% 
% Example 2
% This example matches electrode positions for subject p10 to the wang
% full probability atlas and benson atlas
%     projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual_ecog_recoded'; 
%     subject           = 'p10';
%     atlasName         = {'wang15_fplbl', 'benson14_varea'};
% electrode_table = bidsEcogMatchElectrodesToAtlas(projectDir, subject, [], atlasName);
%
% See also, bidsSpecifySessions, bidsEcogReadElectrodeFile,
%           bidsEcogReadSurfFile, bidsEcogReadAtlasFile, bidsEcogGetMatchedAtlas,
%           bidsEcogPlotElectrodesOnMesh

% IG, BAIR 2020; K.Yuasa, 2022

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
if ~iscell(atlasName), atlasName = {atlasName}; end
atlasName(ismember(atlasName,'none')) = [];

% <optionals>
if ~exist('printSummary','var') && exist('surfaceType','var') && (islogical(surfaceType)||isnumeric(surfaceType))
    printSummary = surfaceType;
    surfaceType  = [];
end
if ~exist('surfaceType','var') && exist('thresh','var') && (ischar(thresh)||isstring(thresh))
    surfaceType  = thresh;
    thresh       = [];
end

% <thresh>
if ~exist('thresh', 'var') || isempty(thresh), thresh = inf; end

% <surfaceType>
if ~exist('surfaceType','var') || isempty(surfaceType)
    % surfaceType should be 'pial' to match electrode locations in BAIR project
    surfaceType = 'pial';
end

% <printSummary>
if ~exist('printSummary', 'var') || isempty(printSummary), printSummary = true; end

%% Read in files and matched ECoG electrodes and MRI atlas

[~,~,~,~,atlasName] = interpretAtlasNames(atlasName);
[matched_atlas_vals, electrode_table, matched_vertices, keep_idx, indices] =...
    bidsEcogGetMatchedAtlas(projectDir, subject, session, atlasName, thresh, surfaceType);

%% Make output

% Add matched nodes to electrode table
electrode_table.node_x(keep_idx) = matched_vertices(:,1);
electrode_table.node_y(keep_idx) = matched_vertices(:,2);
electrode_table.node_z(keep_idx) = matched_vertices(:,3);
electrode_table.node_idx(keep_idx) = indices';

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
            electrode_table.(atlasName{a}) = repmat({'none'}, height(electrode_table),1);
            matched_areas = cell(size(vals));
            idx = vals>0;
            matched_areas(idx) = area_labels(vals(idx));
            matched_areas(~idx) = {'none'};
            electrode_table.(atlasName{a})(keep_idx) = matched_areas;
        else % otherwise, add separate column for each area with atlas vals
            ncols = size(vals,2);
            for ii = 1:ncols
                colName = sprintf('%s_%s', atlasName{a}, area_labels{ii});
                electrode_table.(colName) = zeros(height(electrode_table),1);
                electrode_table.(colName)(keep_idx) = vals(:,ii);
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
