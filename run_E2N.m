% close all;
% check paths
if isempty(which('mrVista')), tbUse('vistasoft'); end
addpath('/Volumes/server/Projects/BAIR/ECoG/code');

%% run script for a specific patient (with plotmesh on) to look at electrodes on brain

% e.g.
pID = 661;
thresh = []; % default is 20 mm
plotmesh = 'yes'; % plot meshes with atlases for each subject: yes or no
plotlabel = 'yes'; % plot electrode labels on mesh: yes or no

retinotopicElecs = electrode_to_nearest_node(pID, thresh, plotmesh, plotlabel);

%% loop matching script across patients (plotmesh off to prevent many figures):

pIDs = [439 496 503 507 516 523 559 560 561 562 563 ...
        564 567 568 569 578 584 590 591 595 596 598 ...
        605 607 608 609 610 615 617 625 626 628 630 ...
        632 633 634 637 638 639 640 645 646 647 648 ...
        651 652 653]; % list of patient IDs to loop over
thresh = []; % default is 20 mm

% run script for each patient
out_all = cell(length(pIDs),1);
for i = 1:length(pIDs)
    out_all{i} = electrode_to_nearest_node(pIDs(i), thresh);
end

%% calculate distribution or electrodes across atlas regions:
atlasNames = {'wang2015_atlas','template_areas'};

patient_count = 0;
for a = 1:length(atlasNames)
    area_names = out_all{1}.(atlasNames{a}).area_names; % first patient in loop should have data
    count_total = zeros(1,length(area_names));
    for i = 1:length(out_all)
        if ~isempty(out_all{i})
            count = out_all{i}.(atlasNames{a}).area_count;
            count_total = count_total+count;
            area_names = out_all{i}.(atlasNames{a}).area_names;
            if a == 1
                patient_count = patient_count+1;
            end
        end
    end
    figure('Name', atlasNames{a});hold on;
    bar(1:length(count_total), count_total, 'k');
    set(gca, 'Ylim', [0 max(count_total)+1]);
    set(gca, 'XLim', [0 length(area_names)+1], 'XTick', 1:1:length(area_names), 'XTickLabel', area_names); 
    xlabel('area');
    ylabel('number of electrodes')
    title(['total number of patients: ' num2str(patient_count)]);
end

disp(['patients with data from list: ' num2str(patient_count) ' out of ' num2str(length(pIDs))]); 



