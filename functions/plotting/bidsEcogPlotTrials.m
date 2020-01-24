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
%                        - channels:  Cell array with channel names to be plotted.
%                                           default: all channels                            
%                        - conditions:Name of stimulus conditions to plot
%                                           default: all conditions                            
%                        - type: single trial, average 
%                        - add_max:
%                        - add_atlas:
%                        - ....
%                        - ....
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
if ~exist('specs', 'var') || isempty(specs)
    specs = struct();
end
if ~isfield(specs,'epoch_t') || isempty(specs.epoch_t)
    specs.epoch_t = [-0.2 1];
end
if ~isfield(specs,'nbase_t') || isempty(specs.epoch_t)
    specs.base_t = [];
end
if ~isfield(specs,'chan_names') || isempty(specs.chan_names)
    specs.chan_names = [];
end
if ~isfield(specs,'stim_names') || isempty(specs.stim_names)
    specs.stim_names = [];
end

% <plot save>
if ~exist('savePlot', 'var') || isempty(savePlot)
    savePlot = true;
end

%% PREPARE DATA

% Load data
    
dataPath = fullfile(projectDir, 'derivatives', inputFolder);
writePath = fullfile(projectDir, 'derivatives', 'ECoGFigures');
                        
[data, channels, events] = bidsEcogGetPreprocData(dataPath, subject, sessions, tasks, runnums, description);
                                          
% Select channels
if isempty(specs.chan_names) 
    chan_idx = contains(channels.type, {'ecog', 'seeg'}); 
else
    chan_idx = ecog_matchChannels(specs.chan_names, channels.name);
end

data = data(chan_idx,:);
channels = channels(chan_idx,:);

% Select trials
if isempty(specs.stim_names)
    specs.stim_names = unique(events.trial_type);
end
stim_inx = cell(length(specs.stim_names, 1));

for ii = 1:length(specs.stim_names)
    if ~isnumeric(specs.stim_names)
        stim_inx{ii} = find(contains(events.trial_name, specs.stim_names{ii}));
    else
        stim_inx{ii} = find(events.trial_name == specs.stim_names(ii));
    end
end

events = events(stim_inx{:},:);

% Epoch the data
[epochs, t] = ecog_makeEpochs(data, events.onset, specs.epoch_t, channels.sampling_frequency(1));  
fprintf('[%s] Step 4: Found %d epochs across %d runs and %d sessions \n', ...
    mfilename, size(epochs,2), length(unique(events.run_name)), length(unique(events.session_name)));

% Baseline correct the epochs
run_idx = [];
switch description
    case 'reref'
        baseType = 'subtractwithintrial';
    case 'broadband'
        baseType = 'percentsignalchange';
end
       
[epochs] = ecog_normalizeEpochs(epochs, t, specs.base_t, baseType, idx);
            
%% MAKE PLOT
                

%% 

if savePlot

   figSaveDir = fullfile(writePath, sprintf('sub-%s', subject), sprintf('ses-%s', session), 'figures');
   if ~exist(figSaveDir, 'dir')
        mkdir(figSaveDir); fprintf('[%s]: Creating a figure directory for sub-%s, ses-%s\n', mfilename, subject, session); 
   end    

   t = ((0:hdr.nSamples-1)/hdr.Fs); 

   fprintf('[%s] Saving Broadband figures to %s \n',mfilename, figSaveDir);

   % Plot the first channel by default
   chan_index = find(contains(lower(channels.type), 'ecog') & contains(channels.status, 'good'));
   channel_plot = chan_index(1);

   figure('Name', sprintf('broadband %s', channels.name{channel_plot}));

   % Plot the time courses before and after CAR
   subplot(2,1,1); hold on
   plot(t,data(channel_plot,:),'k')
   %legend('raw data');
   xlabel('Time (s)'); ylabel('Voltage');set(gca, 'FontSize', 16);
   title(sprintf('%s raw %s',channels.name{channel_plot}, description));

   % Plot the spectra before and after CAR
   subplot(2,1,2),hold on
   plot(t,broadband(channel_plot,:),'k')
   %legend('broadband time course')
   xlabel('Time (s)'); ylabel('Broadband estimate');set(gca, 'FontSize', 16);
   title(sprintf('%s %s \n %s',channels.name{channel_plot}, 'broadband timecourse', methodstr));

   % Generate a name for the figure
   figureName = fullfile(figSaveDir,sprintf('%s_%s',fname_out, channels.name{channel_plot}));
   saveas(gcf, figureName, 'png');

   close all;
end
clear data hdr;

   