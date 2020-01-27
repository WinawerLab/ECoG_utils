
run_idx = {'01','02', '03', '04', '05', '06'};
clear specs;
%for ii = 1:length(run_idx)
    projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual'; 
    subject           = 'som748';
    session           = 'nyuecog01';
    task              = 'temporalpattern';
    runnums           =  [];%run_idx(ii); %, '02', '03', '04'};
    inputFolder       = 'ECoGBroadband';
    description       = 'broadband';
    
    specs.epoch_t     = [-0.2 1.5];
    %specs.chan_names  = {'GB017', 'GB049', 'GB065', 'GB018', 'GB050', 'GB066',  'GB088',  'GB104',  'GB120'};
    specs.chan_names  = {'GB034', 'GB050'};
    %specs.chan_names  = 'GB';
    specs.plot_type   = 'averageSE';
    %specs.stim_names  = {'ONEPULSE-1', 'ONEPULSE-2', 'ONEPULSE-3', 'ONEPULSE-4', 'ONEPULSE-5', 'ONEPULSE-6'};
    specs.stim_names  = {'TWOPULSE-1', 'TWOPULSE-2', 'TWOPULSE-3', 'TWOPULSE-4', 'TWOPULSE-5', 'TWOPULSE-6'};
    bidsEcogPlotTrials(projectDir, subject, session, task, runnums, ...
        inputFolder, description, specs);
%end    