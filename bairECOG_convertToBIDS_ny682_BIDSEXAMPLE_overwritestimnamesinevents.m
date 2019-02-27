%% %% OVERWRITE STIMULUS NAMES in events.tsv

BIDSDataDir = '/Users/winawerlab/matlab/git/bids-examples/';
projectName = 'ieeg_visual_multimodal';
patientID   = 682;
sub_label   = ['som' num2str(patientID)]; 

ses_label = {'nyu3t01', 'somecog01'};
ses_type = {'func', 'ieeg'};

stimDir = '/Volumes/server/Projects/BAIR/Data/BIDS/visual/stimuli';

%% OVERWRITE STIM NAMES IN EVENTS
for ii = 1:length(ses_label)    
    dataWriteDir = fullfile(BIDSDataDir, projectName, sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label{ii}), ses_type{ii});
    EventFiles = dir(fullfile(dataWriteDir,'*events.tsv'));
    for jj = 1:length(EventFiles)
        T = readtable(fullfile(EventFiles(jj).folder,EventFiles(jj).name), 'FileType', 'text');
        % load original stimulus file
        load(fullfile(stimDir, T.stim_file{1}));
        
        % generate list of stimulus IDs
        stimIDs = [];
        for aa = 1:size(stimulus.im_cell,2)
            stimIDsToAdd = 1:size(stimulus.im_cell{aa},3);
            if contains(T.stim_file{1}, {'run-2', 'run-02', 'run-4', 'run-04'})
                stimIDsToAdd = stimIDsToAdd+size(stimulus.im_cell{aa},3);
            end
            stimIDs = [stimIDs stimIDsToAdd];
        end
        
        stim_file_NEW = [];
        for kk = 1:height(T)
            %stimCatNumber = stimulus.cat(kk);
            %stim_inx = stimCatNumber - min(stimulus.cat)+1;
            %stimCatName = stimulus.categories{stim_inx};
            stimCatNumber = num2str(T.trial_type(kk));
            stimCatName = num2str(T.trial_name{kk});
            stimID = stimIDs(T.stim_file_index(kk));
            if stimID < 10
                stimIDName = ['0' num2str(stimID)];
            else
                stimIDName = num2str(stimID);
            end
            stim_file_NEW{kk} = ['stim_' stimCatNumber '_' stimCatName '_' stimIDName '.png'];
        end
      
         %OLD  % overwrite stim_file column
%         stim_file_NEW = [];
%         for kk  = 1:height(T)
%             stim_file_NEW{kk} = ['stim_' num2str(T.trial_type(kk)) '.png'];
%         end
        T.stim_file = stim_file_NEW';
        
        % remove stim_file_index column
        if ii == 1
            T = T(:,1:6);
        elseif ii == 2
            T = T(:,[1:6 8]);
        end
        
        % overwrite file
        writetable(T, fullfile(EventFiles(jj).folder,EventFiles(jj).name), 'FileType','text', 'Delimiter', '\t')
    end
end

%% FIX MISSING ISIs IN ECOG TSVs
ses_label = {'nyu3t01', 'somecog01'};
ses_type = {'func', 'ieeg'};

ii = 1; %fMRI

dataWriteDir = fullfile(BIDSDataDir, projectName, sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label{ii}), ses_type{ii});
EventFiles1 = dir(fullfile(dataWriteDir,'*temporalpattern*events.tsv'));

ii = 2;
dataWriteDir = fullfile(BIDSDataDir, projectName, sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label{ii}), ses_type{ii});
EventFiles2 = dir(fullfile(dataWriteDir,'*temporalpattern*events.tsv'));

for jj = 1:length(EventFiles1)
    T1 = readtable(fullfile(EventFiles1(jj).folder,EventFiles1(jj).name), 'FileType', 'text');    
    T2 = readtable(fullfile(EventFiles2(jj).folder,EventFiles2(jj).name), 'FileType', 'text');    
    T2.ISI = T1.ISI;
    writetable(T2, fullfile(EventFiles2(jj).folder,EventFiles2(jj).name), 'FileType','text', 'Delimiter', '\t')
end



