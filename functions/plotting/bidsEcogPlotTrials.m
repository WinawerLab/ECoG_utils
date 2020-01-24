function bidsEcogPlotTrials(projectDir, subject, sessions, tasks, runnums, ...
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
%                        - epoch_t: Epoch window, default: [-0.5 1];
%                        - base_t:  Baseline time window, e.g. [-0.2 0]. 
%                                           default: all t<0 in epochTime.
%                        - chan_names:  Cell array with channel names to be plotted.
%                                           default: all channels                            
%                        - stim_names:Name of stimulus conditions to plot
%                                           default: all conditions                            
%                        - plot_type: single trial, average, averageSE
%                        - plot_colors:
%                        - plot_ylim: 
%                        - add_max:
%                        - add_atlas:
%                       
%     savePlot:         Flag indicating whether to generate plots of the 
%                       broadband timecourse in derivatives/ECoGfigures
%                           default: true
%
% Example 1
% This example computes broadband timecourses for all the data for all
% sessions, tasks, and runs found for this subject, using default settings 
%     projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual'; 
%     subject           = 'som648';
%     bidsEcogBroadband(projectDir, subject)
% 
% See also bidsSpecifySessions.m bidsSpecifyData.m ecog_extractBroadband.m

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
if ~isfield(specs,'epoch_t') || isempty(specs.epoch_t), specs.epoch_t = [-0.2 1];end
if ~isfield(specs,'base_t') || isempty(specs.base_t), specs.base_t = [min(specs.epoch_t) 0];end
if ~isfield(specs,'chan_names'), specs.chan_names = []; end
if ~isfield(specs,'stim_names'), specs.stim_names = []; end
if ~isfield(specs,'plot_type') || isempty(specs.plot_type), specs.plot_type = 'average'; end
if ~isfield(specs,'plot_cmap') || isempty(specs.plot_cmap), specs.plot_cmap = 'copper'; end
if ~isfield(specs,'plot_ylim'), specs.plot_ylim = []; end

% <plot save>
if ~exist('savePlot', 'var') || isempty(savePlot), savePlot = true; end

%% PREPARE DATA

% Load data
    
dataPath = fullfile(projectDir, 'derivatives', inputFolder);
writePath = fullfile(projectDir, 'derivatives', 'ECoGFigures');
                        
[data, channels, events] = bidsEcogGetPreprocData(dataPath, subject, sessions, tasks, runnums, description);
                                          
% Select channels
if isempty(specs.chan_names) 
    chan_idx = contains(channels.type, {'ecog', 'seeg'}); 
else
    %chan_idx = ecog_matchChannels(specs.chan_names, channels.name);
    chan_idx = contains(channels.name, specs.chan_names);
end

data = data(chan_idx,:);
channels = channels(chan_idx,:);

% Select trials
if isempty(specs.stim_names)
    specs.stim_names = unique(events.trial_type);
end
stim_idx = cell(length(specs.stim_names),1);

for ii = 1:length(specs.stim_names)
    if ~isnumeric(specs.stim_names)
        stim_idx{ii} = find(contains(events.trial_name, specs.stim_names{ii}));
    else
        stim_idx{ii} = find(events.trial_type == specs.stim_names(ii));
    end
end

events = events(vertcat(stim_idx{:}),:);

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
            
%% MAKE PLOT

% Determine how many figures to make
elec_groups = unique(channels.group);

% Get electrode indices for each group; check if this patient has a HD
% grid, if so, reorder channels and make separate plots for bottom and top
group_names = []; 
chan_idx = []; 
c = 1;
for ii = length(elec_groups)
    if contains(elec_groups{ii}, {'gridB', 'HDgrid'})
        chan_idx_HD = find(contains(channels.group, elec_groups{ii}));
        [grid_idx, gridPlotNames] = ecog_sortElecsOnHDGrid(channels.name(chan_idx_HD));
        chan_idx{c} = chan_idx_HD(grid_idx{1});
        chan_idx{c+1} = chan_idx_HD(grid_idx{2});
        group_names{c} = sprintf('%s %s', elec_groups{ii}, gridPlotNames{1});
        group_names{c+1} = sprintf('%s %s', elec_groups{ii}, gridPlotNames{2});
        c = c+2;
    else 
        chan_idx{c}   = find(contains(channels.group, elec_group));
        chan_names{c} = channels.name(chan_idx{c});
        group_names{c} = elec_groups{ii};
    end
end



% Set plot settings
if ischar(specs.plot_cmap)
    cmap = eval(specs.plot_cmap);
    colors = cmap(1:round((length(cmap)/length(stim_idx))):length(cmap),:);
else 
    colors = specs.plot_cmap(1:length(stim_idx),:);
end

nFig = length(group_names);

for ii = 1:nFig
    
    figureName = sprintf('%s %s %s %s', subject, group_names{ii}, tasks, specs.stim_names);
    figure('Name', figureName);
    set(gcf, 'Position', [400 200 1800 1200]);
    
    
    % Determine number of subplots
    % Decide how many subplots are needed
    nPlot = length(chan_idx);
    nRow  = ceil(sqrt(nPlot));
    nCol  = ceil(sqrt(nPlot));
    if nPlot <= (nRow*nCol)-nCol
        nRow = nRow-1;
    end
    
    % Loop over electrodes
    for ee = 1:length(chan_idx)
        
        % Add axis labels and 
        if ee == 1
            legend(specs.stim_names);
            yLabel = sprintf('%s (%s)', description, channels.units{1}); 
        else
            yLabel = [];
        end

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
                    llim = this_trial - std(this_trial,0,2, 'omitnan')/sqrt(length(stim_idx{ss}));
                    ulim = this_trial + std(this_trial,0,2, 'omitnan')/sqrt(length(stim_idx{ss}));
                    CI = [llim ulim];     
            end
            
            % Plot
            ecog_plotSingleTimeCourse(t, this_trial, CI, colors(ss,:), plotTitle, yLabel, specs.plot_ylim, [])
            
        end    
    end
    
    %% save Plot?

    if savePlot

       figSaveDir = fullfile(writePath, sprintf('sub-%s', subject), 'figures');
       if ~exist(figSaveDir, 'dir')
            mkdir(figSaveDir); fprintf('[%s]: Creating a figure directory for sub-%s\n', mfilename, subject); 
       end    

       fprintf('[%s] Saving figures to %s \n',mfilename, figSaveDir);
       saveas(gcf, figureName, 'png');

    end
end




   