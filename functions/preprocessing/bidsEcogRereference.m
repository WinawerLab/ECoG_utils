function bidsEcogRereference(projectDir, subject, sessions, tasks, runnums, outputFolder, savePlot)
% Applies a regression-based common average reference (CAR) to
% bids-formatted ECoG data, and writes out the rereferenced data to an
% output folder in the bids derivatives folder. 
%
% Notes:
% This code uses the 'group' and 'status' columns in channels.tsv to
% determine which channels to average into a common reference. For a given
% group, only 'good' channels are included in the reference average, which
% is then regressed out of all channels of that group. If there is no group
% column, channels are grouped by 'type' (a bids-required column) instead.
% If there is no status column, all channels are assumed to be good.
%
% TO DO: Also write out a projection matrix file? (see draft bids derivative)
% TO DO: If following BIDS spec, outputFolder should be 'pipelinename'
% TO DO: Also write out a README file to explain what this is?
%
% Input
%     projectDir:       path where the BIDS projects lies (string)
%     subject:          BIDS subject name (string, all lower case)
%     sessions:         BIDS session name (string, all lower case)
%                           default: all sessions with 'ecog' in the name
%     tasks:            one or more BIDS tasks (string or cell array of strings)
%                           default: all tasks in session
%     runnums:          BIDS run numbers (vector or cell array of vectors)
%                           default: all runs for specified tasks
%     outputFolder:     Name of folder where rereferenced data is placed
%                           default: 'ECOGpreprocessedCAR'
%     savePlot:         generate plots of timecourse and spectra for a
%                       sample electrode before and after CAR in a separate 
%                       'figures' folder in the derivaties folder 
%                           default: true
%
% Example 1
% This example rereferences all the data for all sessions, tasks, and runs
% found for this subject, and generates plots
%     projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual'; 
%     subject           = 'som726';
%     bidsEcogRereference(projectDir, subject)
% 
% Example 2
% This example rereferences the raw data for all runs of a specific session 
% and tasks, and does not generate any plots
%     projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual'; 
%     subject           = 'som726';
%     session           = 'nyuecog03';
%     task              = 'prf'; 
%     bidsEcogRereference(projectDir, subject, session, task, [], [], 0);


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

% <outputFolder>
if ~exist('outputFolder', 'var') || isempty(outputFolder)
    outputFolder = 'ECOGpreprocessedCAR';
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
    fprintf('[%s] Starting CAR for subject: %s session: %s\n', mfilename, subject, session);
    
    % <writeDir>
    writeDir = fullfile(projectDir, 'derivatives', outputFolder, subject, session);
    if ~exist(writeDir, 'dir')
        mkdir(writeDir); fprintf('[%s] Creating a CAR output folder for sub-%s, ses-%s\n', mfilename, subject, session); 
    end    
    
    sessionDir = fullfile(projectDir, sprintf('sub-%s', subject), sprintf('ses-%s', session));
    if ~exist(sessionDir, 'dir')
        error('data directory not found: %s', sessionDir); 
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
           
           % Read in the json file:   
           jsonReadFile = fullfile(sessionDir, 'ieeg', sprintf('%s_ieeg.json', fname));
           if ~exist(jsonReadFile, 'file'), error('data file not found: %s', dataReadFile); end
           fprintf('[%s] Reading in ieeg json file: %s\n', mfilename, jsonReadFile); 
           ieeg_json = jsonread(jsonReadFile);
           
           % Read in the channels file
           chanFile = fullfile(sessionDir, 'ieeg', sprintf('%s_channels.tsv', fname));
           if ~exist(chanFile, 'file'), error('channels file not found: %s', chanFile); end
           fprintf('[%s] Reading in channels file: %s\n', mfilename, chanFile); 
           channels   = readtable(chanFile, 'FileType', 'text');
           
           % Apply CAR
           [data_reref, channels_reref, group_indices, group_names] = ecog_performCAR(data, channels);           
          
           % Save out the rereferenced data to the derivatives folder
           dataWriteFile = fullfile(writeDir, 'ieeg', sprintf('%s_desc-reref_ieeg.eeg', fname));
           fprintf('[%s] Writing new data file: %s\n', mfilename, dataWriteFile); 
           ft_write_data(dataWriteFile, data_reref, 'header', hdr, 'dataformat', 'brainvision_eeg');
           
           % In the ieeg.json file, update the iEEGreference field with a
           % description of the common average procedure:
           ieeg_json.iEEGReference = 'A common average reference (CAR) was computed separately for each group of good channels and regressed out of each channel.';
           jsonWriteFile = fullfile(writeDir, 'ieeg', sprintf('%s_desc-reref_ieeg.json', fname));
           fprintf('[%s] Writing new ieeg json file: %s\n', mfilename, jsonWriteFile); 
           json_options.indent = '    '; 
           jsonwrite(jsonWriteFile,ieeg_json, json_options)
           
           % Save out the rereferenced channels file to the derivatives folder
           chanWriteFile = fullfile(writeDir, 'ieeg', sprintf('%s_desc-reref_channels.tsv', fname));
           fprintf('[%s] Writing new channels file: %s\n', mfilename, chanWriteFile); 
           writetable(channels_reref,chanWriteFile,'FileType','text','Delimiter','\t');
           
           % Copy over the events.tsv file
           eventsReadFile = fullfile(sessionDir, 'ieeg', sprintf('%s_events.tsv', fname));
           eventsWriteFile = fullfile(writeDir, 'ieeg', sprintf('%s_desc-reref_events.tsv', fname));
           copyfile(eventsReadFile, eventsWriteFile)
           fprintf('[%s] Copying over events file: %s\n', mfilename, eventsWriteFile); 

%% DIAGNOSTICS: Look at the effect of CAR
           if savePlot
               
               figSaveDir = fullfile(writeDir, 'figures');
               if ~exist(figSaveDir, 'dir')
                    mkdir(figSaveDir); fprintf('[%s]: Creating a derivatives figures directory for sub-%s, ses-%s\n', mfilename, subject, session); 
               end    
    
               t = ((0:hdr.nSamples-1)/hdr.Fs); 

               fprintf('[%s] Saving CAR figures to %s \n',mfilename, figSaveDir);

               for ee = 1:length(group_indices)

                   chan_index = group_indices{ee};

                   if ~isempty(chan_index)
                        
                       if isfield(summary(channels), 'status') 
                           good_channels = find(contains(channels(chan_index,:).status, 'good'));
                           channel_plot = chan_index(good_channels(1));
                       else
                           channel_plot = chan_index(1);
                       end
                       
                       figure('Name', sprintf('car regress %s', group_names{ee}));

                       % Plot the time courses before and after CAR
                       subplot(1,2,1); hold on
                       plot(t,data(channel_plot,:),'k')
                       plot(t,data_reref(channel_plot,:),'g')
                       legend({'before CAR','after CAR'})
                       xlabel('Time (s)'); ylabel('Voltage');set(gca, 'FontSize', 16);
                       title(channels.name(channel_plot));

                       % Plot the spectra before and after CAR
                       subplot(1,2,2),hold on
                       [pxx,freqs] = pwelch(data(channel_plot,:)',hdr.Fs,0,hdr.Fs,hdr.Fs);
                       [pxx2,~] = pwelch(data_reref(channel_plot,:)',hdr.Fs,0,hdr.Fs,hdr.Fs);
                       plot(freqs,pxx,'k')
                       plot(freqs,pxx2,'g'); 
                       set(gca, 'YScale', 'log')
                       xlabel('Frequency (Hz)'); ylabel('Log power');set(gca, 'FontSize', 16);
                       title(channels.name(channel_plot));
                        
                       % Generate a name for the figure
                       figureName = fullfile(figSaveDir,sprintf('%s_desc-reref_%s',fname, group_names{ee}));
                       saveas(gcf, figureName, 'png');
                   end
               end
               close all;
           end
           clear data hdr;
       end
    end
    if cleartasks, tasks = []; end 
    if clearrunnums, runnums = []; end 
end
fprintf('[%s] Done with CAR! \n', mfilename); 
