
run_idx = {'01','02', '03', '04', '05', '06'};

for ii = 1:length(run_idx)
    projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual'; 
    subject           = 'som748';
    session           = 'nyuecog04';
    task              = [];
    runnums           = run_idx(ii); %, '02', '03', '04'};
    inputFolder       = 'ECoGCAR';
    description       = 'reref';
    specs.chan_names  = {'GB017', 'GB049', 'GB065', 'GB018', 'GB050', 'GB066',  'GB088',  'GB104',  'GB120'};
    specs.plot_type   = 'averageSE';
    specs.stim_names  = 'faces';

    bidsEcogPlotTrials(projectDir, subject, session, task, runnums, ...
        inputFolder, description, specs);
end
    