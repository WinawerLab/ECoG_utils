
%% NYU patients - line noise at 60 Hz
projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual'; 
subject           = 'som648';
bands             = [[70 80]; [80 90]; [90 100]; [100 110]; [130 140]; [140 150]; [150 160]; [160 170]; [190 200]];
% avoiding harmonics of 60 because of line noise

bidsEcogBroadband(projectDir, subject, [], [], [], [], bands);
     
 
%% UMCU patients - line noise at 50 hz

projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual'; 
subject           = 'beilen';
bands             = [[60 70]; [70 80]; [80 90]; [110 120]; [120 130]; [130 140]; [160 170]; [170 180]; [180 190]];
% avoiding harmonics of 60 because of line noise

bidsEcogBroadband(projectDir, subject, [], [], [], [], bands);
