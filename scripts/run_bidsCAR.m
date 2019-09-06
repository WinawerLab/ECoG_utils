

%% Example 1
% This example rereferences all the data for all sessions, tasks, and runs
% found for this subject, and generates plots
     projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual'; 
     subject           = 'chaam';
     bidsEcogRereference(projectDir, subject)

     
%% Example 2
% This example rereferences the raw data for all runs of a specific session 
% and tasks, and does not generate any plots
     projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual'; 
     subject           = 'som726';
     session           = 'nyuecog03';
     task              = 'prf'; 
     bidsEcogRereference(projectDir, subject, session, task, [], [], 0);
