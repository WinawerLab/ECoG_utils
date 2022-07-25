function [area_labels, area_cmap, units, area_match] = getAtlasLabels(atlasName, isnorm, issmry)

% [area_labels, area_cmap, units, area_match] = GETATLASLABELS(atlasName)
%   Provides labels for a given atlas and provides a colormap
% 
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
%       

% IG, BAIR 2020; K.Yuasa, 2022

% Check flags
if nargin < 1, error('atlasName not defined.'); end
if ~exist('isnorm','var') || isempty(isnorm),  isnorm = false;  end
if ~exist('issmry','var') || isempty(issmry),  issmry = false;  end

% Check atlasName
[atlasName, issmry] = checkAtlasName(atlasName, isnorm, issmry);

% Load atlas labels
switch atlasName
    
    case {'wang2015_atlas','wang15_mplbl','wangarea'} % Wang maximum probability map

            if issmry
            [area_labels, area_cmap, area_match] = wang15_smry_atlas;
            else
            [area_labels, area_cmap] = wang15_atlas;
            end
            units = 'area';
                       
	case {'wang15_fplbl','wangprob'} % Wang full probability map

            if issmry
            [area_labels, ~, area_match] = wang15_smry_atlas;
            else
            area_labels = wang15_atlas;
            end
            area_cmap   = autumn(64);
            units = 'probability';

	case {'benson14_varea','bensonarea'} % Noahs anatomically derived template

            if issmry
            [area_labels, area_cmap, area_match] = benson14_smry_atlas;
            else
            [area_labels, area_cmap] = benson14_atlas;
            end
            units = 'area';

    case {'template_areas','benson20_mplbl','hcparea'} % old Noah templates & Noah's HCP template

            [area_labels, area_cmap] = benson20_atlas;
            units = 'area';
        
    case {'benson20_fplbl','hcpprob'} % Noah's HCP full probability map
        
            area_labels = benson20_atlas;
            area_cmap   = autumn(64);
            units = 'probability';

	case 'glasser16_atlas' % Glasser atlas
   
            [area_labels, area_cmap] = glasser16_label; 
            units = 'parcellation';
            
    otherwise
        
            area_labels = [];
            area_cmap = [];
            units = [];
      
end

if ~exist('area_match','var')
    area_match = 1:length(area_labels);
end
end

%%%%%%%%%%%% Sub function %%%%%%%%%%%%%%
function  [atlasName, issmry] = checkAtlasName(atlasName, isnorm, issmry)
% Correct 'atlasName' in the suitable manner for GETATLASLABELS

% Convert _fplbl_norm to _mplbl
if isnorm
    atlasName = regexprep(atlasName,'_fplbl.*','_mplbl');
else
    atlasName = regexprep(atlasName,'_fplbl_norm.*','_mplbl');
end

% Check smry
issmry = issmry || startsWith(atlasName, 'smry_');
atlasName = strrep(atlasName,'smry_','');
end

