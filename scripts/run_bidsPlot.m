

projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual'; 
subject           = 'som748';
session           = 'nyuecog04';
task              = [];
runnums           = [];
inputFolder       = 'ECoGBroadband';
description       = 'broadband';
specs.chan_names  = 'GB';

bidsEcogPlotTrials(projectDir, subject, session, task, runnums, ...
    inputFolder, description, specs);
