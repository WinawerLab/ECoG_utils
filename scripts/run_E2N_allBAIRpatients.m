

% close all;
if isempty(which('electrode_to_nearest_node'))
    tbUse('ECoG_utils');
end

%% run script for a specific patient (with plotmesh on) to look at electrodes on brain

%pIDs  = [648 645 661 668 674 682 692 708 718 723]; % BAIR patients

% contrasttemporal selection
sub_label = {'648'; 
             '661'; 
             '674'; 
             '692'; 
             '708'; 
             '718';
             '723';
             'umcuchaam';
             'umcubeilen'};
         
% contrasttemporal selection
% sub_label = {'648'; 
%              '661'; 
%              '674'; 
%              '692'; 
%              '708'; 
%              '718';
%              '723'};
         
%contrasttemporal selection
%sub_label = {'umcuchaam'};%
             %'umcubeilen'};

% run script for each patient
specs = [];
out_all = cell(length(sub_label),1);
for i = 1:length(sub_label)
    
    specs.pID           = sub_label{i};
    specs.atlasNames    = {'benson14_varea'};% Can plot only one at a time
    if i == 1
        specs.plotmesh  = 'right';
        figure;hold on;
    else
        specs.plotmesh  = 'none';
    end
    specs.plotelecs     = 'right';
    specs.plotlabel     = 'no';  
    
    out_all{i} = electrode_to_nearest_node_fsaverage2(specs);
end
material dull;
% % run script for each patient
% specs = [];
% out_all = cell(length(pIDs),1);
% for i = 1:length(pIDs)
%     
%     specs.pID           = num2str(pIDs(i));
%     specs.atlasNames    = {'wang15_mplbl'};% Can plot only one at a time
%     specs.plotmesh      = 'both';
%     specs.plotelecs     = 'yes';
%     specs.plotlabel     = 'no';  
%     specs.patientPool   = 'SOM';
%     %specs.elecFile      = sprintf('/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/preprocessed/sub-som%d/som%d_fsaverage_electrode_coordinates.tsv', pIDs(i), pIDs(i));
%     %specs.elecFile      = sprintf('/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/preprocessed/sub-som%d/som%d_fsaverage_electrode_coordinates.tsv', pIDs(i), pIDs(i));
%     specs.thresh = [];
%     
%     out_all{i} = electrode_to_nearest_node_fsaverage(specs);
% end