% test

loadDir = '/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/preprocessed/';
sub_label = {'som648', 'som723'};
ses_label = {'nyuecog01', 'nyuecog01'};

for ii = 1:length(sub_label)
    filename = fullfile(loadDir, sprintf('sub-%s', sub_label{ii}), sprintf('ses-%s', ses_label{ii}), sprintf('sub-%s_ses-%s_epoched.mat', sub_label{ii}, ses_label{ii}));
    load(filename);
    if ii == 1 
        events = trials.events;
    else
        events = [events; trials.events];
    end
end

%load('/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/preprocessed/sub-som648/ses-nyuecog01/sub-som648_ses-nyuecog01_epoched.mat')