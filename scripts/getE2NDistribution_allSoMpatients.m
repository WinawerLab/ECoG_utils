%% loop matching script across patients (plotmesh off to prevent many figures):

% list of patient IDs to loop over
pIDs  = [648 645 661 668 674 682 692 694]; % BAIR patients

pIDs2 = [439 496 503 507 516 523 559 560 561 562 563 ...
        564 567 568 569 578 584 590 591 595 596 598 ...
        605 607 608 609 610 615 617 625 626 628 630 ...
        632 633 634 637 638 639 640 646 647 651 652 ...
        653 660 662 663 664 665 671 673 679]; % SOM patients


pIDs = [pIDs pIDs2];

% run script for each patient
specs = [];
out_all = cell(length(pIDs),1);
for i = 1:length(pIDs)
    
    specs.pID         = num2str(pIDs(i));
    specs.plotmesh  = 'none'; % plot meshes with atlases for each subject: yes or no
  
    out_all{i} = electrode_to_nearest_node(specs);
    %out_all{i} = electrode_to_nearest_node_fsaverage(specs);
end

%% calculate distribution or electrodes across atlas regions:
atlasNames = {'wang2015_atlas','template_areas', 'benson14_varea'};

patient_count = zeros(length(atlasNames),1);
for a = 1:length(atlasNames)
    area_names = out_all{1}.(atlasNames{a}).area_names; % first patient in loop should have data
    count_total = zeros(1,length(area_names));
    for i = 1:length(out_all)
        if ~isempty(out_all{i}) && ~isempty(out_all{i}.(atlasNames{a}))
            count = out_all{i}.(atlasNames{a}).area_count;
            count_total = count_total+count;
            area_names = out_all{i}.(atlasNames{a}).area_names;
            patient_count(a,1) = patient_count(a,1)+1;
        end
    end
    figure('Name', atlasNames{a});hold on;
    bar(1:length(count_total), count_total, 'k');
    set(gca, 'Ylim', [0 max(count_total)+1]);
    set(gca, 'XLim', [0 length(area_names)+1], 'XTick', 1:1:length(area_names), 'XTickLabel', area_names); 
    xlabel('area');
    ylabel('number of electrodes')
    title(['total number of patients: ' num2str(patient_count(a,1))]);
    disp(['patients with electrodes in ' atlasNames{a} ':' num2str(patient_count(a,1)) ' out of ' num2str(length(pIDs))]); 
end


%% identify which patients had any visual coverage

VisCov.names = [];
VisCov.indices = [];

for a = 1:length(atlasNames)
    VisCov.names.(atlasNames{a}) = [];
	VisCov.indices.(atlasNames{a}) = [];
    for p = 1:length(out_all)
        if ~isempty(out_all{p})
            if ~isempty(out_all{p}.(atlasNames{a}))
                if ~isempty(out_all{p}.(atlasNames{a}).elec_labels)
                    VisCov.names.(atlasNames{a}) = [VisCov.names.(atlasNames{a}) {out_all{p}.patientID}];
                    VisCov.indices.(atlasNames{a}) = [VisCov.indices.(atlasNames{a}) p];
                end
            end
        end
    end
end

% To get a list with all the subjects with coverage in e.g. Wang atlas:
test = out_all(VisCov.indices.wang2015_atlas);


