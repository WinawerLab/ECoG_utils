function bidsEcogRereference(projectDir, subject, sessions, tasks, runnums, outputFolder, savePlot)
% Applies a regression-based common average reference (CAR) to
% bids-formatted ECoG data, and writes out the rereferenced data to an
% output folder in the bids derivatives folder. 
%
% Notes
% This code uses the 'group' and 'status' columns in channels.tsv to
% determine which channels to average into a common reference. For a given
% group, only 'good' channels are included in the reference average, which
% is then regressed out of all channels of that group. If there is no group
% column, channels are grouped by 'type' (a bids-required column) instead.
% If there is no status column, all channels are assumed to be good.
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
%     outputFolder:     Name of folder where rereferenced data is placed
%                           default: 'ECoGCAR'
%     savePlot:         generate plots of timecourse and spectra for a
%                       sample electrode before and after CAR in a separate 
%                       'figures' folder in the derivatives folder 
%                           default: true
%
% Example 1
% This example rereferences all the data for all sessions, tasks, and runs
% found for this subject, and generates plots
%     projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual'; 
%     subject           = 'p01';
%     bidsEcogRereference(projectDir, subject)
% 
% Example 2
% This example rereferences the raw data for all runs of a specific session 
% and tasks, and does not generate any plots
%     projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual'; 
%     subject           = 'p03';
%     session           = 'nyuecog01';
%     task              = 'prf'; 
%     bidsEcogRereference(projectDir, subject, session, task, [], [], 0);
%
% See also bidsSpecifySessions.m bidsSpecifyData.m ecog_performCAR.m
%
% IG, BAIR 2019

% <projectDir>
if ~exist('projectDir', 'var') || isempty(projectDir), error('projectDir not defined'); end    

% <subject>
if ~exist('subject', 'var') || isempty(subject), error('subject not defined'); end

% <session>
if ~exist('sessions', 'var') || isempty(sessions)
    [sessions] = bidsSpecifySessions(projectDir, subject);
    idx = find(contains(lower(sessions), {'ecog', 'iemu'}));
    if isempty(idx)
        error('no ECOG sessions found for this subject');
    else
        sessions = sessions(idx);
    end
end

% <outputFolder>
if ~exist('outputFolder', 'var') || isempty(outputFolder), outputFolder = 'ECoGCAR'; end

% <plots>
if ~exist('savePlot', 'var') || isempty(savePlot), savePlot = true; end

%% Check formats and initialize

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
    
    % define paths
    writePath = fullfile(projectDir, 'derivatives', outputFolder);
      
    for jj = 1:length(tasks)
       for kk = 1:length(runnums{jj})
           
           task = tasks{jj};
           runnum = runnums{jj}{kk};
           fprintf('[%s] Task = %s, Run = %s \n', mfilename, task, runnum);
                     
           [data, channels, events, ieeg_json, hdr] = bidsEcogReadFiles(projectDir, subject, session, task, runnum);
           
           % Overwrite channels sampling frequency with header info
           channels.sampling_frequency = repmat(hdr.Fs, [height(channels) 1]);
           
           % Apply CAR
           [data_reref, channels_reref, group_indices, group_names] = ecog_performCAR(data, channels);           
          
           % In the ieeg.json file, update the iEEGreference field
           ieeg_json.iEEGReference = 'A common average reference (CAR) was computed separately for each group of good channels and regressed out of each channel.';
          
           % Update the description and save out the data 
           [fname_out] = bidsEcogWriteFiles(writePath, subject, session, task, runnum, 'reref', ...
                data_reref, channels_reref, events, ieeg_json, hdr);         

%% DIAGNOSTICS: Look at the effect of CAR
           if savePlot
               
               figSaveDir = fullfile(writePath, sprintf('sub-%s', subject), sprintf('ses-%s', session), 'figures');

               if ~exist(figSaveDir, 'dir')
                    mkdir(figSaveDir); fprintf('[%s]: Creating a figure directory for sub-%s, ses-%s\n', mfilename, subject, session); 
               end    
    
               t = ((0:hdr.nSamples-1)/hdr.Fs); 

               fprintf('[%s] Saving CAR figures to %s \n',mfilename, figSaveDir);

               for ee = 1:length(group_indices)

                   chan_index = group_indices{ee};

                   if ~isempty(chan_index)
                        
                       if isfield(summary(channels), 'status') 
                           good_channels = find(contains(channels_reref(chan_index,:).status, 'good'));
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
                       % title(channels.name(channel_plot));

                       % Plot the spectra before and after CAR
                       subplot(1,2,2),hold on
                       nonnanidx = ~isnan(data(channel_plot,:));
                       [pxx,freqs] = pwelch(data(channel_plot,nonnanidx)',hdr.Fs,0,hdr.Fs,hdr.Fs);
                       [pxx2,~] = pwelch(data_reref(channel_plot,nonnanidx)',hdr.Fs,0,hdr.Fs,hdr.Fs);
                       plot(freqs,pxx,'k')
                       plot(freqs,pxx2,'g'); 
                       set(gca, 'YScale', 'log')
                       set(gca, 'XLim', [0 256]) 
                       % maximum freq resolution for typical sample rate at
                       % NYU SOM (UMCU data have higher sample rates);
                       % hardcoded rather than made read from hdr so it's
                       % easier to compare across NYU and UMCU
                       xlabel('Frequency (Hz)'); ylabel('Log power');set(gca, 'FontSize', 16);
                       % title(channels.name(channel_plot));
                        
                       % Generate a name for the figure
                       figureName = fullfile(figSaveDir,sprintf('%s_desc-reref_%s',fname_out, group_names{ee}));
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
