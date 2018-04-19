%% loop matching script across patients (plotmesh off to prevent many figures):

pIDs = [439 496 503 507 516 523 559 560 561 562 563 ...
        564 567 568 569 578 584 590 591 595 596 598 ...
        605 607 608 609 610 615 617 625 626 628 630 ...
        632 633 634 637 638 639 640 646 647  ...
        651 652 653 660 662 663 664 665]; % list of patient IDs to loop over
thresh = []; % default is 20 mm

pIDs2 = [645 648 661];

% run script for each patient
out_all = cell(length(pIDs),1);
for i = 1:length(pIDs)
    
    specs.pID         = num2str(pIDs(i));
    specs.patientPool = 'SOM';
    specs.elecFile    = [];
    specs.fsDir       = ['/Volumes/server/Freesurfer_subjects/som' specs.pID];
    specs.thresh      = []; % default is 20 mm
    specs.plotmesh  = 'no'; % plot meshes with atlases for each subject: yes or no
    specs.plotlabel = 'no'; % plot electrode labels on mesh: yes or no

    out_all{i} = electrode_to_nearest_node(specs);
    %out_all{i} = electrode_to_nearest_node_fsaverage(specs);
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

%% identify which patients had any visual coverage

VisCov.wang      = [];
VisCov.templates = [];
for p = 1:length(out_all)
    if ~isempty(out_all{p})
        if ~isempty(out_all{p}.template_areas.elec_labels)
            VisCov.templates = [VisCov.templates {out_all{p}.patientID}];
        end
        if ~isempty(out_all{p}.wang2015_atlas.elec_labels)
            VisCov.wang = [VisCov.wang {out_all{p}.patientID}];
        end  
    end
end

