function [area_labels, area_cmap] = wang15_atlas(idx)
% WANG15_LABEL outputs the atlas labels and respective colormap of
% wang15_atlas
% 
% Labels come from https://hub.docker.com/r/nben/occipital_atlas/
% 
% Wang L, et al. (2015) Probabilistic Maps of Visual Topography in Human Cortex, Cerebral Cortex.
% https://doi.org/10.1093/cercor/bhu277

% K.Yuasa, BAIR 2022

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

if exist('idx','var') && ~isempty(idx)
    area_labels = area_labels(idx);
    area_cmap   = area_cmap(idx,:);
end
