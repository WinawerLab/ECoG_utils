%% projectDir        
projectDir = '/Volumes/server/Projects/BAIR/Data/BIDS/visual'; 
%%
subjects = {'som648'; 
            'som661'; 
            'som674'; 
            'som682'; 
            'som692'; 
            'som708'; 
            'som709'; 
            'som718'; 
            'som723'; 
            'som726';
            'beilen'; 
            'chaam'};
%subjects = {'beilen'; 
%            'chaam'};
        
for ii = 1:length(subjects)
    subject = subjects{ii};
    bidsEcogRereference(projectDir, subject);
end

%%
subject = 'som748';
session = 'nyuecog04';
bidsEcogRereference(projectDir, subject, session);

%% 
projectDir = '/Volumes/server/Projects/BAIR/Data/BIDS/motor'; 
subject = 'som756';
session = 'nyuecog01';
bidsEcogRereference(projectDir, subject, session);
