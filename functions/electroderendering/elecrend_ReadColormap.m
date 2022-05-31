function [area_cmap, units, col_range] = elecrend_ReadColormap(atlasName)
% [area_cmap, units, col_range] = elecrend_ReadColormap(atlasName)
%   reads colormap data from colormap_atlasName.mat

% K.Yuasa, BAIR 2022

% Set default values for outputs
area_cmap = [];
units     = '';
col_range = [];

% Colormap name
colormap_file = ['colormap_' atlasName '.mat'];

% Read colormap
if exist(colormap_file,'file') && ismember({'cmap'},who('-file',colormap_file))
    
    fprintf('[%s] Loading %s \n',mfilename, colormap_file);
    cmap = load(colormap_file);  cmap = cmap.cmap;
    
    switch atlasName
        case {'benson14_eccen'}
            area_cmap = cmap(15:end-200,2:4);
            col_range = [0 20];
        case {'benson14_sigma'}
            area_cmap = cmap(:,2:4);
            col_range = [0 10];
        case {'benson14_angle'}
            area_cmap = cmap(:,2:4);
            col_range = [0 360];
        otherwise
            area_cmap = cmap(:,2:4);
            col_range = minmax(cmap(:,1)');
    end
    
    if endsWith(atlasName,{'_eccen','_angle','_sigma'})
        units = 'degrees';
    else
        units = 'unknown';
    end
    
elseif ~ismember({'none'}, atlasName)
    error('Colormap for %s is unknown',atlasName);
    
end

end