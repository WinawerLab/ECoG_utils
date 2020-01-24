

projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual'; 
subject           = 'som748';
session           = 'nyuecog01';
task              = 'temporalpattern';
runnums           = [];
inputFolder       = 'ECoGBroadband';
description       = 'broadband';


bidsEcogPlotTrials(projectDir, subject, session, task, runnums, ...
    inputFolder, description);
