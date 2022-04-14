function [area_labels, area_cmap] = benson20_atlas(idx)
% BENSON20_LABEL outputs the atlas labels and respective colormap of
% benson20_atlas
% 
% Benson NC, et al. (2020) Variability of the Surface Area of the V1, V2, and V3 Maps in a Large Sample of Human Observers.
% https://www.biorxiv.org/content/10.1101/2020.12.30.424856v2

% K.Yuasa, BAIR 2022
           
area_labels = {'V1', 'V2', 'V3'}; 

area_cmap   = [255 255 0; 0 255 255; 0 0 255]./255;

if exist('idx','var') && ~isempty(idx)
    area_labels = area_labels(idx);
    area_cmap   = area_cmap(idx,:);
end
