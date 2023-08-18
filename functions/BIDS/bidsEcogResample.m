function bidsEcogResample(projectDir, subject, sessions, tasks, runnums, srate, inputFolder, outputFolder, description)
% Resamples to bids-formatted ECoG data, and writes out the
% resampled data to an output folder in the bids derivatives folder.
% 
% bidsEcogResample(projectDir, subject, [sessions], [tasks], [runnums], ...
%    srate, [inputFolder], [outputFolder], [description])
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
%     srate:            Integer indicating a desired sampling rate in Hz.
%     inputFolder:      Name of a derivatives folder where broadband is
%                           to be computed on (e.g., rereferenced data)
%                           default: ECoGCAR                    
%     outputFolder:     Name of folder where rereferenced data is placed
%                           default: ECoGResampled
%     description:      String stating the 'desc-' label in the name of
%                       the input data files
%                           default: 'reref' 
%
% Example
% Resample all data to 512 Hz for all sessions, tasks, and runs of the subject
%     projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual'; 
%     subject           = 'p01';
%     srate             = 512;
%     bidsEcogRereference(projectDir, subject, [], [], [], srate)
%
% See also bidsSpecifySessions.m bidsSpecifyData.m
%
% KY, BAIR 2023

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

% <srate>
if ~exist('srate', 'var') 
    srate = [];
end
srate = round(srate);

% <inputFolder>
if ~exist('inputFolder', 'var') || isempty(inputFolder), inputFolder = 'ECoGCAR'; end

% <outputFolder>
if ~exist('outputFolder', 'var') || isempty(outputFolder), outputFolder = 'ECoGResampled'; end

% <description>
if ~exist('description', 'var') || isempty(description),   description = 'reref';  end

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
    fprintf('[%s] Resampling initiated at %d Hz for subject: %s, session: %s\n', mfilename, srate, subject, session);

    % define paths
    dataPath = fullfile(projectDir, 'derivatives', inputFolder);
    writePath = fullfile(projectDir, 'derivatives', outputFolder);
      
    for jj = 1:length(tasks)
       for kk = 1:length(runnums{jj})
           
           task = tasks{jj};
           runnum = runnums{jj}{kk};
           fprintf('[%s] Task = %s, Run = %s \n', mfilename, task, runnum);
                     
           [data, channels, events, ieeg_json, hdr] = bidsEcogReadFiles(dataPath, subject, session, task, runnum, description);
           
           if isempty(srate) || srate == hdr.Fs
             warning('[%s] Skipping resampling: Sampling rate not specified or equal to the original rate.', mfilename);
           else
             % Apply resampling
             data = resample(data', srate, hdr.Fs)';

             % In the ieeg.json file, update the iEEGreference field
             ieeg_json.iEEGReference = sprintf('Data resampled to %d Hz from original %d Hz.', srate, hdr.Fs);

             % Overwrite headder and channels sampling frequency
             hdr.Fs = srate;
             channels.sampling_frequency(:) = hdr.Fs;

             % Update the description and save out the data
             [fname_out] = bidsEcogWriteFiles(writePath, subject, session, task, runnum, 'reref', ...
                 data, channels, events, ieeg_json, hdr);
           end
       end
    end
    if cleartasks, tasks = []; end 
    if clearrunnums, runnums = []; end 
end
fprintf('[%s] Done with resampling! \n', mfilename); 
