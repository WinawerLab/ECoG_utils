function  bidsEcogPlotTrials(projectDir, subject, sessions, tasks, runnums, ...
    inputFolder, description, specs, savePlot)

% Plots a time-varying plot of bids-formatted ECoG data; will epoch into
% separate trials, and perform a baseline correction.
%
% bidsEcogPlotTrials(projectDir, subject, [sessions], [tasks], [runnums], ...
%    [inputFolder], [description], [plotSpecs], [savePlot])
%
% Input
%     projectDir:       path where the BIDS project lies (string)
%     subject:          BIDS subject name (string, all lower case)
%     sessions:         BIDS session name (string, all lower case)
%                           default: all sessions with 'ecog' in the name
%     tasks:            one or more BIDS tasks (string or cell array of strings)
%                           default: all tasks in session
%     runnums:          BIDS run numbers (vector or cell array of vectors)
%                           default: all runs for specified tasks
%     inputFolder:      Name of a derivatives folder where broadband is
%                           to be computed on (e.g., rereferenced data)
%                           default: ECoGCAR  
%     description:      String stating the 'desc-' label in the name of
%                       the input data files
%                           default: ...
%     specs:            A struct specifying instructions for what to plot,
%                       with the following fields:
%                        - epoch_t: Epoch window, default: [-0.2 1.2];
%                        - base_t:  Baseline time window, e.g. [-0.2 0]. 
%                                           default: all t<0 in epochTime.
%                        - chan_names: Cell array with channel names to be plotted.
%                                           default: all channels                            
%                        - stim_names: Name of stimulus conditions to plot
%                                           default: all conditions                            
%                        - plot_type: single trial, average, averageSE
%                        - plot_colors: colormap name or RGB array
%                        - plot_ylim: scale of y axis
%                        - add_max:
%                        - add_atlas:
%                       
%     savePlot:         Flag indicating whether to save the plots in 
%                       derivatives/ECoGfigures
%                           default: true
%
% Example
% projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual'; 
% subject           = 'som748';
% session           = 'nyuecog04';
% task              = [];
% runnums           = '01';
% inputFolder       = 'ECoGBroadband';
% description       = 'broadband';
% specs.chan_names  = 'GB';
% 
% bidsEcogPlotTrials(projectDir, subject, session, task, runnums, ...
%     inputFolder, description, specs);
% 
% See also bidsEcogGetPreprocData.m bidsECoGBroadband.m

% <projectDir>
if ~exist('projectDir', 'var') || isempty(projectDir)
    error('projectDir not defined');
end    

% <subject>
if ~exist('subject', 'var') || isempty(subject)
    error('subject not defined');
end

if ~exist('sessions', 'var'), sessions = []; end
if ~exist('tasks', 'var'), tasks = []; end
if ~exist('runnums', 'var'), runnums = []; end

% <inputFolder>
if ~exist('inputFolder', 'var') || isempty(inputFolder)
    inputFolder = 'ECoGBroadband';
end

% <description>
if ~exist('description', 'var') || isempty(description)
    description = 'broadband';
end

% <specs>
if ~exist('specs', 'var') || isempty(specs), specs = struct(); end
if ~isfield(specs,'epoch_t') || isempty(specs.epoch_t), specs.epoch_t = [-0.2 1.2];end
if ~isfield(specs,'base_t') || isempty(specs.base_t), specs.base_t = [min(specs.epoch_t) 0];end
if ~isfield(specs,'chan_names'), specs.chan_names = []; end
if ~isfield(specs,'stim_names'), specs.stim_names = []; end
if ~isfield(specs,'plot_type') || isempty(specs.plot_type), specs.plot_type = 'averageSE'; end
if ~isfield(specs,'plot_cmap') || isempty(specs.plot_cmap), specs.plot_cmap = 'parula'; end
if ~isfield(specs,'plot_ylim'), specs.plot_ylim = []; end
if ~isfield(specs, 'average_stims'), specs.average_stims = 0; end
% <plot save>
if ~exist('savePlot', 'var') || isempty(savePlot), savePlot = true; end

%% PREPARE DATA

% Load data
    
dataPath = fullfile(projectDir, 'derivatives', inputFolder);
writePath = fullfile(projectDir, 'derivatives', 'ECoGFigures');
                        
[data, channels, events] = bidsEcogGetPreprocData(dataPath, subject, sessions, tasks, runnums, description);
                                          
% Select channels
chan_names = specs.chan_names;
if ~iscell(chan_names), chan_names = {chan_names}; end
if isempty([chan_names{:}]) 
    chan_idx = contains(channels.type, {'ecog', 'seeg'}); 
else
    chan_idx = contains(channels.name, chan_names);
end

data = data(chan_idx,:);
channels = channels(chan_idx,:);

% Epoch the data
[epochs, t] = ecog_makeEpochs(data, events.onset, specs.epoch_t, channels.sampling_frequency(1));  
fprintf('[%s] Found %d epochs across %d runs and %d sessions \n', ...
    mfilename, size(epochs,2), length(unique(events.run_name)), length(unique(events.session_name)));

% Baseline correct the epochs
switch description
    case 'reref'
        baseType = 'subtractwithintrial';
    case 'broadband'
        baseType = 'percentsignalchange';
end       
[epochs] = ecog_normalizeEpochs(epochs, t, specs.base_t, baseType);
fprintf('[%s] Baseline correcting epochs using %s \n', mfilename, baseType);
 
% Select trials
stim_names = specs.stim_names;
if ~iscell(stim_names), stim_names = {stim_names}; end
if isempty([stim_names{:}])
    stim_names = unique(events.trial_name);
end
stim_idx = cell(length(stim_names),1);

for ii = 1:length(stim_names)
    if ~isnumeric(stim_names)
        stim_idx{ii} = find(contains(events.trial_name, stim_names{ii}));        
    else
        stim_idx{ii} = find(events.trial_type == stim_names(ii));
    end
end
if specs.average_stims, stim_idx = {vertcat(stim_idx{:})}; stim_names = {[stim_names{:}]}; end

%% MAKE PLOTS

% Set plot settings
if ischar(specs.plot_cmap)
    cmap = eval(specs.plot_cmap);
    colors = cmap(1:round((length(cmap)/length(stim_idx))):length(cmap),:);
else 
    colors = specs.plot_cmap(1:length(stim_idx),:);
end

% Determine how many figures to make
groups = unique(channels.group);
nGroups = length(groups);
chan_groups = []; c = 1;
groupNames = [];

% Get electrode indices for each group 
for ii = 1:nGroups
    
    % Check if group is a HD grid, if so, reorder channels 
    % and make separate plots for bottom and top
    if contains(groups{ii}, {'gridB', 'HDgrid'}) && height(channels) > 64
        grid_idx = find(contains(channels.group, groups{ii}));
        [inx, inxNames] = ecog_sortElecsOnHDGrid(channels,grid_idx);
        chan_groups{c} = inx{1};
        chan_groups{c+1} = inx{2};
        groupNames{c} = sprintf('%s %s', groups{ii}, inxNames{1}); 
        groupNames{c+1} = sprintf('%s %s', groups{ii}, inxNames{2}); 
        c = c+2;
    else      
        chan_groups{c} = find(contains(channels.group, groups{ii}));
        groupNames{c} = groups{ii};
        c = c+1;
    end
end

for ii = 1:length(chan_groups)

    figureName = sprintf('sub %s %s %s %s %s', subject, groupNames{ii}, [stim_names{:}], description, specs.plot_type);

    figure('Name', figureName);
    set(gcf, 'Position', [400 200 1800 1200]);
    
    chan_idx = chan_groups{ii};
    
    % Determine how many subplots to make
    nPlot = length(chan_idx);
    nRow  = ceil(sqrt(nPlot));
    nCol  = ceil(sqrt(nPlot));
    if nPlot <= (nRow*nCol)-nCol
        nRow = nRow-1;
    end
    
    hasLegend = 0;

    % Loop over electrodes
    for ee = 1:length(chan_idx)
        
        if ~isnan(chan_idx(ee))

            subplot(nRow, nCol, ee);
            plotTitle = channels.name(chan_idx(ee));

            % Loop over trial types
            for ss = 1:length(stim_idx)

                this_epoch = epochs(:,stim_idx{ss}, chan_idx(ee));
                CI = [];

                % Determine plot type
                switch specs.plot_type

                    case 'singletrial'
                        this_trial = this_epoch; 

                    case 'average'                    
                        this_trial = mean(this_epoch,2, 'omitnan');

                    case 'averageSE'
                        this_trial = mean(this_epoch,2, 'omitnan');
                        llim = this_trial - std(this_epoch,0,2, 'omitnan')/sqrt(length(stim_idx{ss}));
                        ulim = this_trial + std(this_epoch,0,2, 'omitnan')/sqrt(length(stim_idx{ss}));
                        CI = [llim ulim];     
                end

                % Plot
                ecog_plotSingleTimeCourse(t, this_trial, CI, colors(ss,:), [], [], specs.plot_ylim);
            end
            c = c+1;
            
            % Add axis labels and legend
            if ~hasLegend
                legend(stim_names); hasLegend = 1;
                ylabel(sprintf('%s (%s)', description, channels.units{1}));
                %xlabel('Time (s)');
            end
            
            title(plotTitle);
            
        end
    end
    
    %% save Plot?

    if savePlot

       figSaveDir = fullfile(writePath, sprintf('sub-%s', subject), 'figures');
       if ~exist(figSaveDir, 'dir')
            mkdir(figSaveDir); fprintf('[%s]: Creating a figure directory for sub-%s\n', mfilename, subject); 
       end    

       fprintf('[%s] Saving figures to %s \n',mfilename, figSaveDir);
       saveas(gcf, fullfile(figSaveDir, figureName), 'png');

    end
end




   