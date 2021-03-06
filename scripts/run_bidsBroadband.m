%% projectDir        
projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual'; 

%% NYU patients - line noise at 60 Hz


bands             = [[70 80]; [80 90]; [90 100]; [100 110]; [130 140]; [140 150]; [150 160]; [160 170]; [190 200]];
% avoiding harmonics of 60 because of line noise

subjects = {'som648'; 
            'som661'; 
            'som674'; 
            'som682'; 
            'som692'; 
            'som708'; 
            'som709'; 
            'som718'; 
            'som723'; 
            'som726'};
    
for ii = 1:length(subjects)
    subject = subjects{ii};
    bidsEcogBroadband(projectDir, subject, [], [], [], bands);
end

%function bidsEcogBroadband(projectDir, subject, sessions, tasks, runnums, ...
%    bands, method, inputFolder, outputFolder, description, savePlot)

%%
subject = 'som748';
session = 'nyuecog04';
bidsEcogBroadband(projectDir, subject, session, [], [], bands);

projectDir = '/Volumes/server/Projects/BAIR/Data/BIDS/motor'; 
subject = 'som756';
session = 'nyuecog01';
bidsEcogBroadband(projectDir, subject, session, [], [], bands);

%% UMCU patients - line noise at 50 hz

bands             = [[60 70]; [70 80]; [80 90]; [110 120]; [120 130]; [130 140]; [160 170]; [170 180]; [180 190]];
% avoiding harmonics of 50 because of line noise

subjects = {'beilen'; 
            'chaam'};
        
for ii = 1:length(subjects)
    subject = subjects{ii};
    bidsEcogBroadband(projectDir, subject, [], [], [], bands);
end

     