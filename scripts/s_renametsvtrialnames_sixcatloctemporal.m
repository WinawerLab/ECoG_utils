% s_renametsvtrialnames_sixcatloctemporal

%%% WARNING: WILL OVERWRITE EXISTING EVENT FILES %%%

projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual'; 
subject           = 'som748';
session           = 'nyuecog04';
task              = 'sixcatloctemporal';
overwrite         = 1;

dataDir = fullfile(projectDir, sprintf('sub-%s',subject), sprintf('ses-%s', session), 'ieeg');
d = dir(fullfile(dataDir, sprintf('*_task-%s*events.tsv', task)));


for ii = 1:length(d)
    fname = fullfile(d(ii).folder, d(ii).name);
    events = readtable(fname, 'FileType', 'text');
    
    durs = unique(events.duration);
    
    newtrialnames = [];
    for jj = 1:height(events)
        category_name = upper(events.trial_name{jj});
        this_ISI = events.ISI(jj);
        this_dur = events.duration(jj);
        % determine condition 
        if this_ISI == 0
            condition_name = 'ONEPULSE';
            condition_type = find(this_dur == durs);
        else
            condition_name = 'TWOPULSE';
            condition_type = find(this_ISI == durs);
        end
        newtrialnames{jj} = sprintf('%s-%d-%s', condition_name, condition_type,category_name);
    end
    
    % Overwrite existing tsv file
    events.trial_name = newtrialnames';
    if overwrite
        writetable(events,fname,'FileType','text','Delimiter','\t');
    end

end