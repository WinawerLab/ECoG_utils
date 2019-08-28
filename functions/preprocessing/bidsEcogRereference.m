function bidsEcogRereference(projectDir, subject, sessions, tasks, runnums, makeplot)
% Apply a regression-based common average reference (CAR) to bids-formatted
% ECoG data, and write out the rereferenced data to a folder called
% 'ECOGpreprocessedCAR' in the bids derivaties folder.
%
% Note: This code will use the status column in the channels.tsv file to
% read in which channels are good, and apply CAR only to those channels
% (only good channels are included in the reference, and only those
% channels are rereferenced). If there is no status column, it will be
% assumed that all channels are good and all will be included (in that case
% data should not contain DC, status channels or ECG channels).

% Input
%     projectDir:       path where the BIDS projects lies (string)
%     subject:          BIDS subject name (string, all lower case)
%     sessions:         BIDS session name (string, all lower case)
%                           default: all sessions with 'ecog' in the name
%     tasks:            one or more BIDS tasks (string or cell array of strings)
%                           default: all tasks in session
%     runnums:          BIDS run numbers (vector or cell array of vectors)
%                           default: all runs for specified tasks
%     makeplot:         generate plots of timecourse and spectra for a
%                       sample electrode before and after CAR in a separate 
%                       'figures' folder in the derivaties folder 
%                           default: true
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
%     bidsEcogRereference(projectDir, subject, session, task, [], 0);


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

if ~exist('makeplot', 'var') || isempty(makeplot)
    makeplot = true;
end

if ~iscell(sessions), sessions = {sessions}; end
if ~exist('tasks', 'var'), tasks = []; end
if ~exist('runnums', 'var'), runnums = []; end

%% Perform CAR for each session, tasks and runnums 

for mm = 1:length(sessions)   
    
    [session, tasks, runnums] = bidsSpecifyData(projectDir, subject, sessions{mm}, tasks, runnums);
    
    % <writeDir>
    writeDir = fullfile(projectDir, 'derivatives', 'ECOGpreprocessedCAR', subject, session);
    if ~exist(writeDir, 'dir')
        mkdir(writeDir); fprintf('[%s] Creating a writeDir for sub-%s, ses-%s\n', mfilename, subject, session); 
    end    
    
    sessionDir = fullfile(projectDir, sprintf('sub-%s', subject), sprintf('ses-%s', session));
    if ~exist(sessionDir, 'dir')
        error('data directory not found: %s', sessionDir); 
    end
   
    for jj = 1:length(tasks)
       for kk = 1:length(runnums{jj})
           
           fname = sprintf('sub-%s_ses-%s_task-%s_run-%s', subject, session, tasks{jj}, runnums{jj}{kk});
           chanFile = fullfile(sessionDir, 'ieeg', sprintf('%s_channels.tsv', fname));
           
           % Read in the channels file
           if ~exist(chanFile, 'file')
                error('channels file not found: %s', chanFile); 
           end
           fprintf('[%s] Reading in channels file: %s\n', mfilename, chanFile); 
           channels   = readtable(chanFile, 'FileType', 'text');

           % Read in the data file
           dataReadFile = fullfile(sessionDir, 'ieeg', sprintf('%s_ieeg.eeg', fname));
           if ~exist(dataReadFile, 'file')
                error('data file not found: %s', dataReadFile); 
           end
           fprintf('[%s] Reading in data file: %s\n', mfilename, dataReadFile); 
           hdr = ft_read_header(dataReadFile);
           data = ft_read_data(dataReadFile);
           
           % Apply CAR
           [data_reref, INX, INXNames] = ecog_performCAR(data, channels);           
          
           % Save out the rereferenced data to the derivatives folder
           dataWriteFile = fullfile(writeDir, sprintf('%s_ieeg.eeg', fname));
           fprintf('[%s] Writing new data file: %s\n', mfilename, dataWriteFile); 
           ft_write_data(dataWriteFile, data_reref, 'header', hdr, 'dataformat', 'brainvision_eeg');

           % Save out a log/json file with car description/settings
           
%            % TO DO
%            car_json = [];
%            json_options = [];
%            jsonWriteFile = fullfile(writeDir, sprintf('%s_car.json', fname));    
%            jsonwrite(jsonWriteFile,car_json,json_options);
%            

%% DIAGNOSTICS: Look at the effect of CAR
           if makeplot
               
               figSaveDir = fullfile(writeDir, 'figures');
               if ~exist(figSaveDir, 'dir')
                    mkdir(figSaveDir); fprintf('[%s]: Creating a derivatives figures directory for sub-%s, ses-%s\n', mfilename, subject, session); 
               end    
    
               t = ((0:hdr.nSamples-1)/hdr.Fs); 

               fprintf('[%s] Saving CAR figures...\n',mfilename);

               for mm = 1:length(INX)

                   chan_index = INX{mm};

                   if ~isempty(chan_index)

                        good_channels = find(contains(channels(chan_index,:).status, 'good'));
                        channel_plot = chan_index(good_channels(1));

                        figure('Name', sprintf('car regress %s', INXNames{mm}));

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

                        %set(gcf, 'Position',[1000 700 1100 600]); 
                        
                        % Generate a name for the figure
                        figureName = fullfile(figSaveDir,sprintf('%s_CAR-%s',fname, INXNames{mm}));
                        saveas(gcf, figureName, 'png');
                   end
               end
               close all;
           end
       end
   end
end