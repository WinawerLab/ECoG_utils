function  bidsEcogPlotSpectra(projectDir, subject, sessions, tasks, runnums, ...
    inputFolder, description, specs, savePlot)

% Plots spectral plot of bids-formatted ECoG data; will epoch into
% separate trials.
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
%                                           default: [-1 2];
%                        - fft_stim_t: FFT time window (s)
%                                           default: [0.05 0.55].
%                        - fft_blank_t: FFT time window (s) relative to stimulus
%                                           onset in which NO stimulus was
%                                           presented, to plot as reference
%                                           default: [-0.55 -0.05];
%                        - fft_w: FFT hanning window size (s)
%                                           default: 0.1
%                        - fft_ov: FFT hanning window overlap (s)
%                                           default: 0.05
%                        - reg_erp: Regress out ERP (boolean)
%                                           default: 0.
%                        - chan_names: Cell array with channel names
%                                           default: all channels in channels table
%                        - stim_names: Cell array with stimulus names
%                                           default: all unique names in events table
%                        - plot_type: single trial, average, averageSE
%                        - plot_includeblank: boolean
%                        - plot_colors: colormap name or RGB array
%                        - plot_ylim: limits of y axis
%                                           default: automatic by matlab
%                        - plot_xlim: limits of x axis (frequencies)
%                                           default: [1 200];
%                        - plot_xscale: scale of x axis (string)
%                                           default: 'linear';
%                        - plot_yscale: scale of y axis (string)
%                                           default: 'log';
%                        - average_stims: average all stim_names together
%                                           default: false
%                        - add_atlas: TO DO (read atlas from derivatives/freesurfer)
%                        - fig_subplotdims:
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
    inputFolder = 'ECoGCAR';
end

% <description>
if ~exist('description', 'var') || isempty(description)
    description = 'reref';
end

% <specs>
if ~exist('specs', 'var') || isempty(specs), specs = struct(); end
if ~isfield(specs,'epoch_t') || isempty(specs.epoch_t), specs.epoch_t = [-1 2]; end
if ~isfield(specs,'fft_w') || isempty(specs.fft_w), specs.fft_w = 0.1; end
if ~isfield(specs,'fft_ov') || isempty(specs.fft_ov), specs.fft_ov = 0.05; end
if ~isfield(specs,'fft_stim_t') || isempty(specs.fft_stim_t), specs.fft_stim_t = [0.05 0.55]; end
if ~isfield(specs,'fft_blank_t') || isempty(specs.fft_blank_t), specs.fft_blank_t = [-0.55 -0.05]; end
if ~isfield(specs,'reg_erp') || isempty(specs.reg_erp), specs.reg_erp = 0; end
if ~isfield(specs,'chan_names'), specs.chan_names = []; end
if ~isfield(specs,'stim_names'), specs.stim_names = []; end
if ~isfield(specs,'plot_includeblank'), specs.plot_includeblank = []; end
if ~isfield(specs,'plot_type') || isempty(specs.plot_type), specs.plot_type = 'averageSE'; end
if ~isfield(specs,'plot_cmap') || isempty(specs.plot_cmap), specs.plot_cmap = 'parula'; end
if ~isfield(specs,'plot_xlim') || isempty(specs.plot_xlim), specs.plot_xlim = [1 200]; end
if ~isfield(specs,'plot_ylim'), specs.plot_ylim = []; end
if ~isfield(specs,'plot_xscale') || isempty(specs.plot_xscale), specs.plot_xscale = 'linear';end
if ~isfield(specs,'plot_yscale') || isempty(specs.plot_yscale), specs.plot_yscale = 'log';end

if ~isfield(specs, 'average_stims'), specs.average_stims = 0; end
if ~isfield(specs,'subplotdims') || isempty(specs.subplotdims), specs.subplotdims = []; end
if ~isfield(specs,'subplotidx') || isempty(specs.subplotidx), specs.subplotidx = []; end

% <plot save>
if ~exist('savePlot', 'var') || isempty(savePlot), savePlot = true; end

%% PREPARE DATA

% Load data
dataPath = fullfile(projectDir, 'derivatives', inputFolder);
writePath = fullfile(projectDir, 'derivatives', 'ECoGFigures');

[data, channels, events] = bidsEcogGetPreprocData(dataPath, subject, sessions, tasks, runnums, description);

% Select channels
chan_names = specs.chan_names;
includeChansInFigureName = 0; % if 0, uses the group name only

if ~iscell(chan_names), chan_names = {chan_names}; end
if isempty([chan_names{:}])
    chan_idx = contains(lower(channels.type), {'ecog', 'seeg'});
elseif any(contains(chan_names, string(0:10)))
    % specs.chan_names contains numbers, so the user is (probably)
    % referencing individual channels. Match each individual channel:
    chan_idx = ecog_matchChannels(chan_names, channels.name);
    includeChansInFigureName = 1;
else
    % specs.chan_names does not contain any numbers, so the user is
    % (probably) trying to plot a group of channels based on a common
    % character (e.g. G): match all channels with this character:
    chan_idx = contains(channels.name, chan_names);
    % ALTERNATIVE: use channel_group as second way to select channels
end
if ~any(chan_idx), error('Did not find any matching channels! Please check channel names.'), end

data = data(chan_idx,:);
channels = channels(chan_idx,:);

% Epoch the data
[epochs, t] = ecog_makeEpochs(data, events.onset, specs.epoch_t, channels.sampling_frequency(1));
fprintf('[%s] Found %d epochs across %d runs and %d sessions \n', ...
    mfilename, size(epochs,2), length(unique(events.run_name)), length(unique(events.session_name)));

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
if specs.average_stims, stim_idx = {vertcat(stim_idx{:})}; stim_names = {[stim_names{:}]}; end


%% COMPUTE SPECTRA

% FFT SETTINGS
fft_w = window(@hann,1000*specs.fft_w); % window width (ms)
fft_ov = 1000*specs.fft_ov; % % overlap (ms)
reg_erp = specs.reg_erp; % regress out ERP?

% TRIALS
epochs = permute(epochs, [3 2 1]); % channels x epochs x time

% TIME
stims = ones(1,size(epochs,2));
fft_t = t > specs.fft_stim_t(1) & t < specs.fft_stim_t(2);
srate = channels.sampling_frequency(1);

[f,spectra] = ecog_spectra(epochs,stims,fft_w,fft_t,fft_ov,srate,reg_erp);

hao niif isempty(specs.plot_includeblank) || specs.plot_includeblank
% check inputs
if min(specs.fft_blank_t) < min(specs.epoch_t)
    if isempty(specs.plot_includeblank)
        specs.plot_includeblank = 0;
    else
        error('The FFT time window for blank is smaller than than epoch window. Please adjust or do not plot blank using specs.plot_includeblank = 0.');
    end
else
    specs.plot_includeblank = 1;
    fft_t = t > specs.fft_blank_t(1) & t < specs.fft_blank_t(2);
    [f,spectra_blank] = ecog_spectra(epochs,stims,fft_w,fft_t,fft_ov,srate,reg_erp);

    % add blanks to data and stimulus indices
    nTrials = size(spectra,2);
    spectra = cat(2, spectra, spectra_blank);
    stim_idx = [{(nTrials+1:2*nTrials)'}; stim_idx];
    stim_names = [{'BLANK'} stim_names];
end
end

spectra = permute(spectra,[3 2 1]); % frequencies x epochs x channels

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

if ~iscell(tasks), tasks = {tasks}; end

for ii = 1:length(chan_groups)

    if includeChansInFigureName
        chanFigureNamePart = [chan_names{:}];
    else
        chanFigureNamePart = groupNames{ii};
    end
    figureName = sprintf('sub %s %s %s %s %s %s', subject, chanFigureNamePart, [tasks{:}], [stim_names{:}], 'spectra', specs.plot_type);

    figure('Name', figureName);

    chan_idx = chan_groups{ii};

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
    else
        %set(gcf, 'Position', [400 200 1800 1200/2]);
        set(gcf, 'Position', [0 400 1200 1000]);
    end

    hasLegend = 0;

    % Loop over electrodes
    for ee = 1:length(chan_idx)

        if ~isnan(chan_idx(ee))

            subplot(nRow, nCol, plot_idx(ee));
            plotTitle = channels.name(chan_idx(ee));

            % Loop over trial types
            for ss = 1:length(stim_idx)

                this_spectrum = spectra(:,stim_idx{ss}, chan_idx(ee));
                CI = [];

                % Determine plot type
                switch specs.plot_type

                    case 'singletrial'
                        this_spectrum = this_spectrum;
                        hasLegend = 1; % do not plot legend for single trials

                    case 'average'
                        this_spectrum = mean(this_spectrum,2, 'omitnan');

                    case 'averageSE'
                        this_spectrum = mean(this_spectrum,2, 'omitnan');
                        llim = this_spectrum - std(this_spectrum,0,2, 'omitnan')/sqrt(length(stim_idx{ss}));
                        ulim = this_spectrum + std(this_spectrum,0,2, 'omitnan')/sqrt(length(stim_idx{ss}));
                        CI = [llim ulim];
                end

                % Plot
                ecog_plotSingleSpectrum(f, this_spectrum, CI, colors(ss,:), [], specs.plot_ylim, specs.plot_yscale,[], specs.plot_xlim, specs.plot_xscale);
            end
            c = c+1;

            % Add axis labels and legend
            if ~hasLegend
                legend(stim_names, 'Location', 'best'); hasLegend = 1;
                %ylabel(sprintf('%s (%s)', description, channels.units{1}));
                %xlabel('Time (s)');
            end
            if nPlot <= 4
                ylabel('Spectral power');
                xlabel('Frequency (Hz)');
            end
            setsubplotaxes();
            set(gca, 'FontSize', 20);
            title(plotTitle);
        end
    end

    %% save Plot?

    if savePlot

        figSaveDir = fullfile(writePath, sprintf('sub-%s', subject), 'figures');
        if ~exist(figSaveDir, 'dir')
            mkdir(figSaveDir); fprintf('[%s] Creating a figure directory for sub-%s\n', mfilename, subject);
        end

        fprintf('[%s] Saving figures to %s \n',mfilename, figSaveDir);
        saveas(gcf, fullfile(figSaveDir, figureName), 'png');

    end
end
