
run_idx = {'01','02', '03', '04', '05', '06'};
clear specs;
savePlot = 0;

%for ii = 1:length(run_idx)
    projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual'; 
    subject           = 'som748';
    inputFolder       = 'ECoGBroadband';
    description       = 'broadband';
    
    specs.epoch_t     = [-0.2 1.5];
    %specs.chan_names  = {'GB017', 'GB049', 'GB065', 'GB018', 'GB050', 'GB066',  'GB088',  'GB104',  'GB120'};
    %specs.chan_names  = {'GB034', 'GB050'};
    specs.chan_names  = {'PT0'};
    specs.plot_type   = 'averageSE';
    specs.stim_names  = {'ONEPULSE-1', 'ONEPULSE-2', 'ONEPULSE-3', 'ONEPULSE-4', 'ONEPULSE-5', 'ONEPULSE-6'};
    %specs.stim_names  = {'TWOPULSE-1', 'TWOPULSE-2', 'TWOPULSE-3', 'TWOPULSE-4', 'TWOPULSE-5', 'TWOPULSE-6'};
    %specs.stim_names  = {'FACES', 'HOUSES'};
    specs.plot_ylim  = [-1 5];
    session           = 'nyuecog01';
    task              = 'temporalpattern';%'spatialobject';
    runnums           =  [];%run_idx(ii); %, '02', '03', '04'};
    
    bidsEcogPlotTrials(projectDir, subject, session, task, runnums, ...
        inputFolder, description, specs, savePlot);
    
    %specs.stim_names  = {'TWOPULSE-1-FACES', 'TWOPULSE-2-FACES', 'TWOPULSE-3-FACES', 'TWOPULSE-4-FACES', 'TWOPULSE-5-FACES', 'TWOPULSE-6-FACES'};

    session           = 'nyuecog04';
    task              = [];
    runnums           = [];%run_idx(ii); %, '02', '03', '04'};
    specs.chan_names  = {'GB102', 'GB118', 'PT01', 'PT02'};
    specs.plot_ylim  = [-1 15];

    %specs.stim_names  = {'ONEPULSE-1-FACES', 'ONEPULSE-2-FACES', 'ONEPULSE-3-FACES', 'ONEPULSE-4-FACES', 'ONEPULSE-5-FACES', 'ONEPULSE-6-FACES'};
    specs.stim_names  = {'TWOPULSE-1-FACES', 'TWOPULSE-2-FACES', 'TWOPULSE-3-FACES', 'TWOPULSE-4-FACES', 'TWOPULSE-5-FACES', 'TWOPULSE-6-FACES'};

    %specs.stim_names  = {'FACES', 'BODIES', 'BUILDINGS', 'SCENES'};
    bidsEcogPlotTrials(projectDir, subject, session, task, runnums, ...
        inputFolder, description, specs, savePlot);
% end    