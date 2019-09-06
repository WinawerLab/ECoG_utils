function bidsEcogBroadband(projectDir, subject, sessions, tasks, runnums, ...
    method, bands, inputFolder, outputFolder, savePlot)
% Computes a time-varying broadband measure on for bids-formatted ECoG
% data, and writes out the broadband time course data of equal length to
% the input data to an output folder in the bids derivatives folder.
%
% Notes:
% % TO DO Save out a readme or json that includes information about the bands and
           % methodstr
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
%     method:           method for computing broadband, formatted as a
%                       function handle (see ecog_extractBroadband.m)
%                           default: @(bp) geomean(abs(hilbert(bp)).^2)
%     bands:            band-pass filter frequencies, formatted as a matrix 
%                       (number of bands x 2) or a cell array {[lb,ub], width}, 
%                       (see ecog_extractBroadband.m)
%                           default: {[70 200], 20}
%     inputFolder:      Name of a data folder where broadband is
%                           to be computed on (e.g., rereferenced data)
%                           default: 'ECOGpreprocessedCAR'                            
%     outputFolder:     Name of folder where broadband data is placed
%                           default: 'ECOGpreprocessedBroadband'
%     savePlot:         generate plots of timecourse and spectra for a
%                       sample electrode before and after CAR in a separate 
%                       'figures' folder in the derivatives folder 
%                           default: true
%
% Example 1
% This example rereferences all the data for all sessions, tasks, and runs
% found for this subject, and generates plots
%     projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual'; 
%     subject           = 'som726';
%     bidsEcogBroadband(projectDir, subject)
% 
% Example 2
% This example rereferences the raw data for all runs of a specific session 
% and tasks, and does not generate any plots
%     projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual'; 
%     subject           = 'som726';
%     session           = 'nyuecog03';
%     task              = 'prf'; 
%     bidsEcogBroadband(projectDir, subject, session, task, [], [], 0);
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

% <inputFolder>
if ~exist('inputFolder', 'var') || isempty(outputFolder)
    inputFolder = 'ECOGpreprocessedCAR';
end

% <outputFolder>
if ~exist('outputFolder', 'var') || isempty(outputFolder)
    outputFolder = 'ECOGpreprocessedBroadband';
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

%% Perform CAR for each session, tasks and runnums 

for ii = 1:length(sessions)   
    
    [session, tasks, runnums] = bidsSpecifyData(projectDir, subject, sessions{ii}, tasks, runnums);
    fprintf('[%s] Starting Broadband for subject: %s session: %s\n', mfilename, subject, session);
    
    % <writeDir>
    writeDir = fullfile(projectDir, 'derivatives', outputFolder, subject, session);
    if ~exist(writeDir, 'dir')
        mkdir(writeDir); fprintf('[%s] Creating a Broadband output folder for sub-%s, ses-%s\n', mfilename, subject, session); 
    end    
    
    sessionDir = fullfile(projectDir, 'derivatives', inputFolder, sprintf('sub-%s', subject), sprintf('ses-%s', session));
    if ~exist(sessionDir, 'dir')
        error('input folder not found: %s', sessionDir); 
    end
   
    for jj = 1:length(tasks)
       for kk = 1:length(runnums{jj})
           
           fname = sprintf('sub-%s_ses-%s_task-%s_run-%s', subject, session, tasks{jj}, runnums{jj}{kk});
    
           % Read in the data file
           dataReadFile = fullfile(sessionDir, 'ieeg', sprintf('%s_ieeg.eeg', fname));
           if ~exist(dataReadFile, 'file'), error('data file not found: %s', dataReadFile); end
           fprintf('[%s] Reading in data file: %s\n', mfilename, dataReadFile); 
           hdr = ft_read_header(dataReadFile);
           data = ft_read_data(dataReadFile);
                                 
           % Compute broadband
           [broadband, methodstr] = ecog_extractBroadband(data, hdr.Fs, method, bands);           
          
           % Save out the rereferenced data to the derivatives folder
           dataWriteFile = fullfile(writeDir, 'ieeg', sprintf('%s_desc-broadband_ieeg.eeg', fname));
           fprintf('[%s] Writing new data file: %s\n', mfilename, dataWriteFile); 
           ft_write_data(dataWriteFile, broadband, 'header', hdr, 'dataformat', 'brainvision_eeg');
           
           % In the ieeg.json file, update the iEEGreference field with a
           % description of the common average procedure:
           ieeg_json.iEEGReference = 'A common average reference (CAR) was computed separately for each group of good channels and regressed out of each channel.';
           jsonWriteFile = fullfile(writeDir, 'ieeg', sprintf('%s_desc-broadband_ieeg.json', fname));
           fprintf('[%s] Writing new ieeg json file: %s\n', mfilename, jsonWriteFile); 
           json_options.indent = '    '; 
           jsonwrite(jsonWriteFile,ieeg_json, json_options)
           
           % Copy over the ieeg.json file
           jsonReadFile = fullfile(sessionDir, 'ieeg', sprintf('%s_ieeg.json', fname));
           jsonWriteFile = fullfile(writeDir, 'ieeg', sprintf('%s_desc-broadband_ieeg.json', fname));
           copyfile(jsonReadFile, jsonWriteFile)
           fprintf('[%s] Copying over ieeg json file: %s\n', mfilename, jsonWriteFile); 
           
           % Copy over the channels.tsv file
           chanReadFile = fullfile(sessionDir, 'ieeg', sprintf('%s_channels.tsv', fname));
           chanWriteFile = fullfile(writeDir, 'ieeg', sprintf('%s_desc-broadband_channels.tsv', fname));
           copyfile(chanReadFile, chanWriteFile)
           fprintf('[%s] Copying over channels file: %s\n', mfilename, chanWriteFile); 
           
           % Copy over the events.tsv file
           eventsReadFile = fullfile(sessionDir, 'ieeg', sprintf('%s_events.tsv', fname));
           eventsWriteFile = fullfile(writeDir, 'ieeg', sprintf('%s_desc-broadband_events.tsv', fname));
           copyfile(eventsReadFile, eventsWriteFile)
           fprintf('[%s] Copying over events file: %s\n', mfilename, eventsWriteFile); 
          

%% DIAGNOSTICS: Inspect the broadband time courses
           if savePlot
               
               figSaveDir = fullfile(writeDir, 'figures');
               if ~exist(figSaveDir, 'dir')
                    mkdir(figSaveDir); fprintf('[%s]: Creating a figure directory for sub-%s, ses-%s\n', mfilename, subject, session); 
               end    
    
               t = ((0:hdr.nSamples-1)/hdr.Fs); 

               fprintf('[%s] Saving Broadband figures to %s \n',mfilename, figSaveDir);
                
               % Plot the first channel by default
               channel_plot = 1;

               figure('Name', sprintf('broadband %s', group_names{ee}));

               % Plot the time courses before and after CAR
               subplot(2,1,1); hold on
               plot(t,data(channel_plot,:),'k')
               legend('raw data');
               xlabel('Time (s)'); ylabel('Voltage');set(gca, 'FontSize', 16);
               title(hdr.label(channel_plot));

               % Plot the spectra before and after CAR
               subplot(2,1,2),hold on
               plot(t,broadband(channel_plot,:),'r')
               legend('broadband time course')
               xlabel('Time (s)'); ylabel('Broadband estimate');set(gca, 'FontSize', 16);
               title(hdr.label(channel_plot));
               
%                [pxx,freqs] = pwelch(data(channel_plot,:)',hdr.Fs,0,hdr.Fs,hdr.Fs);
%                [pxx2,~] = pwelch(data_reref(channel_plot,:)',hdr.Fs,0,hdr.Fs,hdr.Fs);
%                plot(freqs,pxx,'k')
%                plot(freqs,pxx2,'g'); 
%                set(gca, 'YScale', 'log')
%                xlabel('Frequency (Hz)'); ylabel('Log power');set(gca, 'FontSize', 16);
%                title(channels.name(channel_plot));

               % Generate a name for the figure
               figureName = fullfile(figSaveDir,sprintf('%s_desc-broadband_%s',fname, group_names{ee}));
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
