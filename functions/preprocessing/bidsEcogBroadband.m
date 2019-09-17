function bidsEcogBroadband(projectDir, subject, sessions, tasks, runnums, ...
    bands, method, inputFolder, outputFolder, description, savePlot)
% Computes a time-varying broadband measure on for bids-formatted ECoG
% data, and writes out the broadband time course data of equal length to
% the input data to an output folder in the bids derivatives folder.
%
% Notes:
% Since one typically wants to do some preprocessing on the raw data (e.g.
% rereferencing) before computing broadband, the inputdata is assumed to
% be located in the derivatives folder and have a 'desc' field in the name.
%
% TO DO Save out a readme?
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
%     bands:            band-pass filter frequencies, formatted as a matrix 
%                       (number of bands x 2) or a cell array {[lb,ub], width}, 
%                       (see ecog_extractBroadband.m)
%                       default: {[60 200], 20}
%     method:           method for computing broadband, formatted as a
%                       function handle (see ecog_extractBroadband.m)
%                           default: @(bp) geomean(abs(hilbert(bp)).^2)
%     inputFolder:      Name of a data folder where broadband is
%                           to be computed on (e.g., rereferenced data)
%                           default: 'ECoGCAR'                            
%     outputFolder:     Name of folder where broadband data is placed
%                           default: 'ECoGBroadband'
%     description:      String stating the 'desc-' label in the name of
%                       the input data files
%                           default: 'reref' 
%     savePlot:         generate plots of timecourse and spectra for a
%                       sample electrode before and after CAR in a separate 
%                       'figures' folder in the derivatives folder 
%                           default: true
%
% Example 1
% This example computes broadband timecourses for all the data for all
% sessions, tasks, and runs found for this subject, using default settings 
%     projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual'; 
%     subject           = 'som648';
%     bidsEcogBroadband(projectDir, subject)
% 
% Example 2
% This example rereferences all the data for all sessions, tasks, and runs
% found for this subject, using custom bands to avoid harmonics of 60Hz
%     projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual'; 
%     subject           = 'som648';
%     bands             = [[70 80]; [80 90]; [90 100]; [100 110]; [130 140]; [140 150]; [150 160]; [160 170]; [190 200]];
%     bidsEcogBroadband(projectDir, subject, [], [], [], bands);
%
% Example 3
% This example rereferences all the data for all sessions, tasks, and runs
% found for this subject, computing broadband amplitude rather than power
%     projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual'; 
%     subject           = 'som648';
%     method            = @(bp)geomean(abs(hilbert(bp)))
%     bidsEcogBroadband(projectDir, subject, [], [], [], [], method);
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

% <session>
if ~exist('sessions', 'var') || isempty(sessions)
    [sessions] = bidsSpecifySessions(projectDir, subject);
    idx = find(contains(lower(sessions), 'ecog'));
    if isempty(idx)
        error('no ECOG sessions found for this subject');
    else
        sessions = sessions(idx);
    end
end

% <method>
if ~exist('method', 'var') 
    method = [];
end

% <bands>
if ~exist('bands', 'var') 
    bands = [];
end

% <inputFolder>
if ~exist('inputFolder', 'var') || isempty(inputFolder)
    inputFolder = 'ECoGCAR';
end

% <outputFolder>
if ~exist('outputFolder', 'var') || isempty(outputFolder)
    outputFolder = 'ECoGBroadband';
end

% <description>
if ~exist('description', 'var') || isempty(description)
    description = 'reref';
end

% <plots>
if ~exist('savePlot', 'var') || isempty(savePlot)
    savePlot = true;
end

if ~iscell(sessions), sessions = {sessions}; end
if ~exist('tasks', 'var'), tasks = []; end
if ~exist('runnums', 'var'), runnums = []; end

cleartasks = 0;
clearrunnums = 0;

if isempty(tasks), cleartasks = 1; end % determines whether tasks will be cleared later in loop
if isempty(runnums), clearrunnums = 1; end

%% Compute broadband for each session, tasks and runnums 

for ii = 1:length(sessions)   
    
    [session, tasks, runnums] = bidsSpecifyData(projectDir, subject, sessions{ii}, tasks, runnums);
    fprintf('[%s] Starting Broadband extraction for sub-%s, ses-%s\n', mfilename, subject, session);
    
    % define paths
    dataPath = fullfile(projectDir, 'derivatives', inputFolder);
    writePath = fullfile(projectDir, 'derivatives', outputFolder);
   
    for jj = 1:length(tasks)
       for kk = 1:length(runnums{jj})
           
           task = tasks{jj};
           runnum = runnums{jj}{kk};
           fprintf('[%s] Task = %s, Run = %s \n', mfilename, task, runnum);
                     
           [data, channels, events, ieeg_json, hdr] = bidsEcogReadFiles(dataPath, subject, session, task, runnum, description);
                                          
           % COMPUTE BROADBAND
           
           % Apply only to those channels that actually have data
           chan_index = find(contains(lower(channels.type), {'ecog', 'seeg'}));
           fprintf('[%s] Found %d ecog and seeg channels, applying broadband to these channels only. \n', mfilename, length(chan_index)); 
           data_bb = data;
           [broadband, methodstr, bandsused] = ecog_extractBroadband(data(chan_index,:), hdr.Fs, method, bands);           
           data_bb(chan_index,:) = broadband;
                       
           % Add columns to channels file about broadband bands and method,
           % update units
           % UNITS (heuristic)
           if contains(methodstr, '.^2')
               unitName = 'uV^2';
           end
           channels.units(chan_index) = {unitName};
           % LOW_CUTOFF 
           if isnumeric(channels.low_cutoff), channels.low_cutoff = num2cell(channels.low_cutoff);end
           if isnumeric(channels.high_cutoff), channels.high_cutoff = num2cell(channels.high_cutoff);end
           channels.low_cutoff(chan_index) = {max(bandsused(:))};
           channels.high_cutoff(chan_index) = {min(bandsused(:))};
           % BANDS, METHODS, BANDWIDTH
           channels.bb_method = repmat({'n/a'}, [height(channels),1]);
           channels.bb_bandwidth = repmat({'n/a'}, [height(channels),1]);
           %channels.bb_bands = repmat({'n/a'}, [height(channels),1]);
           channels.bb_method(chan_index) = {methodstr};
           channels.bb_bandwidth(chan_index) = {diff(bandsused(1,:))};
           %channels.bb_bands(chan_index) = {mat2str(bandsused)};
           
           % Update the description and save out the data 
           [fname_out] = bidsEcogWriteFiles(writePath, subject, session, task, runnum, 'broadband', ...
                data_bb, channels, events, ieeg_json, hdr);         

%% DIAGNOSTICS: Inspect the broadband time courses
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
       end
    end
    if cleartasks, tasks = []; end 
    if clearrunnums, runnums = []; end 
end
fprintf('[%s] Done with Broadband! \n', mfilename); 
