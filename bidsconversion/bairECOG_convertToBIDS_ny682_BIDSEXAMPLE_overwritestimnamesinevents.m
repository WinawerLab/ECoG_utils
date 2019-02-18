%% %% OVERWRITE STIMULUS NAMES in events.tsv

BIDSDataDir = '/Users/winawerlab/matlab/git/bids-examples/';
projectName = 'ieeg_visual_multimodal';
patientID   = 682;
sub_label   = ['som' num2str(patientID)]; 

ses_label = {'nyu3t01', 'somecog01'};
ses_type = {'func', 'ieeg'};

for ii = 1:length(ses_label)    
    dataWriteDir = fullfile(BIDSDataDir, projectName, sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label{ii}), ses_type{ii});
    D = dir(fullfile(dataWriteDir,'*events.tsv'));
    for jj = 1:length(D)
        T = readtable(fullfile(D(jj).folder,D(jj).name), 'FileType', 'text');
        % remove stim_file_index column
        if ii == 1
            T = T(:,1:5);
        elseif ii == 2
            T = T(:,[1:6 8]);
        end
        % overwrite stim_file column
        stim_file_NEW = [];
        for kk  = 1:height(T)
            stim_file_NEW{kk} = ['stim_' num2str(T.trial_type(kk)) '.png'];
        end
        T.stim_file = stim_file_NEW';
        % overwrite file
        writetable(T, fullfile(D(jj).folder,D(jj).name), 'FileType','text', 'Delimiter', '\t')
    end
end
