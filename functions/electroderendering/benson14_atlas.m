function [area_labels, area_cmap] = benson14_atlas(idx)
% BENSON14_LABEL outputs the atlas labels and respective colormap of
% benson14_atlas
% 
% Labels come from Noah email:  values in order from 1-12: V1 V2 V3 hV4 VO1 VO2 LO1 LO2 TO1 TO2 V3b V3a
% See also https://github.com/noahbenson/neuropythy
% 
% Benson NC, et al. (2014) Correction of Distortion in Flattened Representations of the Cortical Surface Allows Prediction of V1-V3 Functional Organization from Anatomy, PLOS Computational Biology. 
% https://doi.org/10.1371/journal.pcbi.1003538

% K.Yuasa, BAIR 2022

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

if exist('idx','var') && ~isempty(idx)
    area_labels = area_labels(idx);
    area_cmap   = area_cmap(idx,:);
end
