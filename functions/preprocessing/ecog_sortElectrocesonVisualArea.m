function elecnames_sorted = ecog_sortElectrodesonVisualArea(elecnames,viselec);

%viselec = trials.viselec.benson14_varea;
%elecnames = whichElectrodes;

sortorder = nan(size(elecnames));

for ii = 1:length(elecnames)
    %inx = contains(viselec.elec_labels,elecnames{ii});
    inx = strcmp(elecnames{ii},viselec.elec_labels);
    if max(inx) > 0
        matchedArea = viselec.area_labels{inx};
        sortorder(ii) = find(strcmp(matchedArea,viselec.area_names));
    end
end

[~,I] = sort(sortorder);
elecnames_sorted = elecnames(I);

end