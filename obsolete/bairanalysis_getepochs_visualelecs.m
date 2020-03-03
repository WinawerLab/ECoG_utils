function [out] = bairanalysis_getepochs_visualelecs(dataLoc, sub_label, ses_label)

% LOAD PREPROCESSED DATA

dataNm  = fullfile(dataLoc, sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label),  sprintf('sub-%s_ses-%s_epoched.mat', sub_label, ses_label));
%a{k}    = load(dataNm);
fprintf('[%s] Loading preprocessed data file %s... \n', mfilename, dataNm);
load(dataNm);

% Collect info on matched electrodes
viselec_wang = []; visarea_wang = [];
for wangname = {'wang2015_atlas', 'wang15_mplbl'}
    if ~isempty(trials.viselec.(wangname{:}))
        viselec_wang = [viselec_wang trials.viselec.(wangname{:}).elec_labels];
        visarea_wang = [visarea_wang trials.viselec.(wangname{:}).area_labels];
    end
end
viselec_benson = []; visarea_benson = [];
if ~isempty(trials.viselec.benson14_varea)
    viselec_benson = [viselec_benson trials.viselec.benson14_varea.elec_labels];
    visarea_benson = [visarea_benson trials.viselec.benson14_varea.area_labels];
end

% Get visual elec indices
[viselec] = unique([viselec_wang viselec_benson])';
name = cell(size(viselec));
wang = repmat({'none'},(size(viselec))); 
% need to put something because empty cells cause trouble with searching
% for specific areas later
benson = repmat({'none'},(size(viselec)));
elInx = []; 
for ii = 1:length(viselec)
    elInx(ii) = ecog_matchChannels(viselec{ii}, trials)';
    name{ii} = viselec{ii};
    wangInx = strmatch(viselec{ii}, viselec_wang, 'exact');
    if ~isempty(wangInx)
        wang{ii} = visarea_wang{wangInx};
    end
    bensonInx = strmatch(viselec{ii}, viselec_benson, 'exact');
    if ~isempty(bensonInx)
        benson{ii} = visarea_benson{bensonInx};
    end
end
visareas = table(name, wang, benson);

% Select visual electrodes
trials.broadband = trials.broadband(elInx,:,:);
trials.evoked = trials.evoked(elInx,:,:);
trials.channels = trials.channels(elInx,:);
blank_trials.broadband = blank_trials.broadband(elInx,:,:);
blank_trials.evoked = blank_trials.evoked(elInx,:,:);
blank_trials.channels = blank_trials.channels(elInx,:);

% Add area matches to channel table
trials.channels.wangatlas = wang;
trials.channels.bensonatlas = benson;
blank_trials.channels.wangatlas = wang;
blank_trials.channels.bensonatlas = benson;

out.sub = sub_label;
out.ses = ses_label;
out.trials = trials;
out.blank_trials = blank_trials;
out.visareas = visareas;

end