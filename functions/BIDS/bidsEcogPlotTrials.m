function  [out] = bidsEcogPlotTrials(projectDir, subject, sessions, tasks, runnums, ...
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
%                       with the following possible fields:
%                        - epoch_t: Epoch window (s)
%                                           default: [-0.2 1.2]
%                        - base_t:  Baseline time window (s)
%                                           default: all t<0 in epochTime.
%                        - chan_names: Cell array with channel names 
%                                           default: all channels in channels table                           
%                        - stim_names: Cell array with stimulus names 
%                                           default: all unique names in events table                         
%                        - plot_type: How to plot the trials (string):
%                                       - 'singletrial': all single trials
%                                       - 'average': trial average per stim
%                                       - 'averageSE': trial average per stim
%                                       with confidence intervals
%                                           default: 'averageSE'                            
%                        - plot_cmap: colormap name (string) or RGB array
%                                           default: 'parula'  
%                        - plot_smooth: extent of smoothing of time course
%                                           default: 0
%                        - plot_ylim: limits of y axis e.g. [-1 10]
%                                           default: automatic by matlab
%                        - plot_includelegend: flag to include legend
%                                           default: true
%                        - average_stims: flag indicating to average across 
%                                         all cells in stim_names (boolean)
%                                           default: false                            
%                        - add_atlas: TO DO (read atlas from derivatives/freesurfer)
%                        - fig_subplotdims: [nrow ncol] array indicating
%                                           layout for the subplots
%                                         all cells in stim_names (boolean)
%                        - fig_subplotidx: 
%                       
%     savePlot:         Flag indicating whether to save the plots in 
%                       derivatives/ECoGfigures
%                           default: true
%
% Example
% projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual_ecog_recoded'; 
% subject           = 'p11';
% session           = 'nyuecog01';
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
if ~isfield(specs,'plot_smooth'), specs.plot_smooth = 0; end
if ~isfield(specs,'plot_includelegend'), specs.plot_includelegend = 1; end
if ~isfield(specs,'average_stims'), specs.average_stims = 0; end
if ~isfield(specs,'subplotdims') || isempty(specs.subplotdims), specs.subplotdims = []; end
if ~isfield(specs,'subplotidx') || isempty(specs.subplotidx), specs.subplotidx = []; end
if ~isfield(specs,'plot_offset') || isempty(specs.plot_offset), specs.plot_offset = 0; end
if ~isfield(specs,'str') || isempty(specs.str), specs.str = []; end

% <plot save>
if ~exist('savePlot', 'var') || isempty(savePlot), savePlot = true; end

%% PREPARE DATA

% Load data
    
dataPath = fullfile(projectDir, 'derivatives', inputFolder);
writePath = fullfile(projectDir, 'derivatives', 'ECoGFigures');
                        
[data, channels, events] = bidsEcogGetPreprocData(dataPath, subject, sessions, tasks, runnums, description, 512);
                                          
% Select channels
chan_names = specs.chan_names;
includeChansInFigureName = 0; % if 0, uses the group name only

if ~iscell(chan_names), chan_names = {chan_names}; end
if isempty([chan_names{:}]) 
    chan_idx = contains(lower(channels.type), {'ecog', 'seeg'}); 
    useSeparateFiguresForGroups = 1;
elseif any(contains(chan_names, string(0:10))) 
    % specs.chan_names contains numbers, so the user is (probably)
    % referencing individual channels. Match each individual channel:
    chan_idx = ecog_matchChannels(chan_names, channels.name);
    includeChansInFigureName = 1;
    useSeparateFiguresForGroups = 0;
else
    % specs.chan_names does not contain any numbers, so the user is
    % (probably) trying to plot a group of channels based on a common
    % character (e.g. G): match all channels with this character:
    % ALTERNATIVE: use channel_group as second way to select channels?
    chan_idx = contains(channels.name, chan_names);
	chan_idx = chan_idx | matches(channels.group, chan_names);
    includeChansInFigureName = 1;
    useSeparateFiguresForGroups = 1;    
end
if ~any(chan_idx), error('Did not find any matching channels! Please check channel names.'), end

data = data(chan_idx,:);
channels = channels(chan_idx,:);

% Epoch the data
if specs.plot_offset % epoch by trial offset

    event_offset = zeros(height(events), 1);

    % Convert trial_name to string if needed
    if iscell(events.trial_name)
        trial_names = string(events.trial_name);
    else
        trial_names = events.trial_name;
    end

    % Identify 'ONE-PULSE' and 'TWO-PULSE' trials
    is_one_pulse = contains(trial_names, "ONE-PULSE");
    is_two_pulse = contains(trial_names, "TWO-PULSE");

    % Compute offsets
    event_offset(is_one_pulse) = events.onsets(is_one_pulse) + events.duration(is_one_pulse);
    event_offset(is_two_pulse) = events.onsets(is_two_pulse) + 2 * events.duration(is_two_pulse) + events.ISI(is_two_pulse);

    % Add to the table
    events.event_offset = event_offset;
    [epochs, t] = ecog_makeEpochs(data, events.event_offset, specs.epoch_t, channels.sampling_frequency(1));
    fprintf('[%s] Found %d epochs across %d runs and %d sessions \n', ...
        mfilename, size(epochs,2), length(unique(events.run_name)), length(unique(events.session_name)));

    % Add offset to figure filename and figure title
    specs.str = 'Relative to offset';

else % epoch by trial onset
    [epochs, t] = ecog_makeEpochs(data, events.onset, specs.epoch_t, channels.sampling_frequency(1));
    fprintf('[%s] Found %d epochs across %d runs and %d sessions \n', ...
        mfilename, size(epochs,2), length(unique(events.run_name)), length(unique(events.session_name)));
end

% Baseline correct the epochs
switch description
    case 'reref'
        baseType = 'subtractwithintrial';
    case 'broadband'
        baseType = 'percentsignalchange';
end
fprintf('[%s] Baseline correcting epochs using %s \n', mfilename, baseType);
[epochs] = ecog_normalizeEpochs(epochs, t, specs.base_t, baseType);
 
% Select trials
fprintf('[%s] Selecting stimulus conditions... \n', mfilename);
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
    fprintf('[%s] Found %d trials for condition %s \n', mfilename, length(stim_idx{ii}), stim_names{ii});
end
if specs.average_stims, stim_idx = {vertcat(stim_idx{:})}; stim_names = {'all stim averaged'}; end%{[stim_names{:}]};

%% MAKE PLOTS

out = [];

% Set plot settings
if ischar(specs.plot_cmap)
    cmap = eval(specs.plot_cmap);
    colors = cmap(1:round((length(cmap)/length(stim_idx))):length(cmap),:);
else 
    colors = specs.plot_cmap(1:length(stim_idx),:);
end

% Determine how many figures to make
chan_groups = []; 
groupNames = [];
if useSeparateFiguresForGroups
    groups = unique(channels.group);
    nGroups = length(groups);
    c = 1;
    % Get electrode indices for each group 
    for ii = 1:nGroups

        % Check if group is a HD grid, if so, reorder channels 
        % and make separate plots for bottom and top
        % TO DO --> replace this by Ken's plotGrid function (no longer split)
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
else
    chan_groups{1} = 1:height(channels);
end

if ~iscell(tasks), tasks = {tasks}; end

for ii = 1:length(chan_groups)
    
    if includeChansInFigureName
        chanFigureNamePart = [chan_names{:}];
    else
        chanFigureNamePart = groupNames{ii};
    end
    figureName = sprintf('sub %s %s %s %s %s %s %s', subject, chanFigureNamePart, [tasks{:}], [stim_names{:}], description, specs.plot_type, specs.str);

    figure('Name', figureName);
    
    chan_idx = chan_groups{ii};
    
    out{ii}.t = t;
    out{ii}.channels = channels(chan_idx(~isnan(chan_idx)),:);
    out{ii}.stims = stim_names;
    out{ii}.colors = colors;
    
    % Determine how many subplots to make and in which order
    nPlot = length(chan_idx);
    if isempty(specs.subplotdims)
        nRow  = ceil(sqrt(nPlot));
        nCol  = ceil(sqrt(nPlot));
        if nPlot <= (nRow*nCol)-nCol, nRow = nRow-1; end
    else
        nRow = specs.subplotdims(1); nCol = specs.subplotdims(2);
    end
    if isempty(specs.subplotidx)
        plot_idx = 1:nPlot;
    else
        plot_idx = specs.subplotidx;
    end
   
    if nPlot > 4
        %set(gcf, 'Position', [400 200 1800 1200]);
        set(gcf, 'Position', get(0, 'Screensize'));
    elseif nPlot < 3
        set(gcf, 'Position', [0 400 1200 500]);
    else
        %set(gcf, 'Position', [400 200 1800 1200/2]);
        set(gcf, 'Position', [0 400 1200 1000]);
    end

    hasLabel = 0;
    hasLegend = 0;
    
    % Loop over electrodes
    for ee = 1:length(chan_idx)
        
        if ~isnan(chan_idx(ee))

            subplot(nRow, nCol, plot_idx(ee));
            plotTitle = channels.name(chan_idx(ee));

            % Loop over trial types
            for ss = 1:length(stim_idx)

                this_epoch = epochs(:,stim_idx{ss}, chan_idx(ee));
                CI = [];

                % Determine plot type
                 switch specs.plot_type

                    case 'singletrial'
                        this_trial = this_epoch; 
                        hasLabel = 1; % do not plot legend for single trials
                        hasLegend = 1;
                        
                    case 'average'                    
                        this_trial = mean(this_epoch,2, 'omitnan');

                    case 'averageSE'
                        this_trial = mean(this_epoch,2, 'omitnan');
                        llim = this_trial - std(this_epoch,0,2, 'omitnan')/sqrt(length(stim_idx{ss}));
                        ulim = this_trial + std(this_epoch,0,2, 'omitnan')/sqrt(length(stim_idx{ss}));
                        CI = [llim ulim]; 
                        %CI = [];
                end

                % Plot
                if specs.plot_smooth > 0
                    this_trial = smooth(this_trial,specs.plot_smooth);
                    if ~isempty(CI)
                        CI(:,1) = smooth(CI(:,1),specs.plot_smooth);
                        CI(:,2) = smooth(CI(:,2),specs.plot_smooth);
                    end
                end
                ecog_plotSingleTimeCourse(t, this_trial, CI, colors(ss,:), [], [], specs.plot_ylim);
                
                % Collect data in output
%                out{ii}.ts(:,ss,ee) = this_trial;
%                out{ii}.ci(:,:,ss,ee) = CI;
            end
            
            % Add axis labels
            if ~hasLabel && specs.plot_includelegend 
                hasLabel = 1;
                % ylabel(sprintf('%s (%s)', description, channels.units{1}));
                ylabel(sprintf('%s (X-fold increase)', description));
                xlabel('Time (s)');
            end

            if nPlot <= 4
                % ylabel(sprintf('%s (%s)', description, channels.units{1}));
                ylabel(sprintf('%s (X-fold increase)', description));
                xlabel('Time (s)');
            end
            setsubplotaxes();
            set(gca, 'FontSize', 20);
            title(plotTitle);          
        end

        if nPlot > 4
            % Create an extra subplot for the legend with invisible lines
            if ee == length(chan_idx)
                subplot(nRow, nCol, plot_idx(ee) + 1);
                hold on

                for ss = 1:length(stim_idx)
                    % Plot but make lines invisible
                    ecog_plotSingleTimeCourse(t, this_trial, CI, colors(ss,:), [], [], specs.plot_ylim);
                end
                set(gca, 'Visible', 'off');
                % Add legend for all trial types
                legend(stim_names, 'Location', 'best');
                axis off; % Hide axes
            end

        elseif nPlot<= 4 && ~hasLegend
            legend(stim_names, 'Location', 'best','FontSize',10);    
            hasLegend = 1;
        end
    end
    
    %% save Plot?

    if savePlot

       figSaveDir = fullfile(writePath, sprintf('sub-%s', subject), 'figures',specs.plot_type);
       if ~exist(figSaveDir, 'dir')
            mkdir(figSaveDir); fprintf('[%s] Creating a figure directory for sub-%s\n', mfilename, subject); 
       end    

       fprintf('[%s] Saving figures to %s \n',mfilename, figSaveDir);
       saveas(gcf, fullfile(figSaveDir, figureName), 'png');

    end
end




   