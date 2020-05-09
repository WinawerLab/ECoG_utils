function [area_labels, area_cmap] = getAtlasLabels(atlasName)

% Provides labels for a given atlas and provides a colormap
% currently included: wang, benson
% to do: glasser

switch atlasName
    
    case {'wang2015_atlas', 'wang15_mplbl'} % Wang maximum probability map

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
            
	case 'benson14_varea' % Noahs anatomically derived template

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

	case 'template_areas' % old Noah templates

            area_labels = {'V1', 'V2', 'V3'}; 
            area_cmap   = [255 255 0; 0 255 255; 0 0 255]./255;

	case 'glasser16_atlas' % Glasser atlas

            area_labels = num2cell(1:180);
            
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
    otherwise
        
        area_labels = [];
        area_cmap = [];
            
end
end

