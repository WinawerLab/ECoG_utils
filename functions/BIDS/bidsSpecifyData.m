function [session, tasks, runnums, modality] = bidsSpecifyData(projectDir, subject,...
    session, tasks, runnums)
% Specify tasks and run numbers and verify paths for BIDS session
% [session, tasks, runnum] = bidsSpecifyData(projectDir, subject, session, [tasks], [runnum])
%
% Input
%     projectDir:       path where the BIDS projects lies (string)
%     subject:          BIDS subject name (string, all lower case)
%     session:          BIDS session name (string, all lower case)
%     tasks:            one or more BIDS tasks (string or cell array of strings)
%                           default: all tasks in session
%     runnums:          BIDS run numbers (vector or cell array of vectors)
%                           default: all runs for specified tasks
%
%
% Output
%     session:          String. See input
%     tasks:            Cell array of strings. See input
%     runnums:          Cell array of vectors. See input 
%     modality:         String indicating data modality. 
%
% Example 1
%     projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual'; 
%     subject           = 'p01';
%     session           = 'nyuecog01';
%     [session, tasks, runnum] = bidsSpecifyData(projectDir, subject, session)
%
% See also bidsSpecifySessions.m


% <projectDir>
if ~exist('projectDir', 'var') || isempty(projectDir)
    error('projectDir not defined');
end    

% <subject>
if ~exist('subject', 'var') || isempty(subject)
    error('subject not defined');
end

% <subjectDir>
subjectDir = fullfile(projectDir, sprintf('sub-%s', subject));
if ~exist(subjectDir, 'dir')
    error('subject dir not found: %s', subjectDir); 
end

% <sessionDir>
sessionDir = fullfile(subjectDir, sprintf('ses-%s', session));
if ~exist(sessionDir, 'dir')
    error('session dir not found: %s', sessionDir); 
end

%% Set the optional inputs

% <tasks>
% MRI and MEG data will have different BIDS formatting, so figure out which
% we're doing. 
if ~exist('tasks', 'var') || isempty(tasks)
    d = dir(fullfile(sessionDir, 'func', '*bold.nii*'));
    modality = 'mri';
    if isempty (d)
        d = dir(fullfile(sessionDir, 'meg', '*meg.sqd'));
        modality = 'meg';
        if isempty (d)
            d = dir(fullfile(sessionDir, 'ieeg', '*ieeg.eeg'));
            modality = 'ecog';
        end
    end
    taskname = cell(1,length(d));
    for ii = 1:length(d)
        taskname{ii} = bidsGet(d(ii).name, 'task');
    end
    tasks = unique(taskname);
end
if ~iscell(tasks), tasks = {tasks}; end

% <runnum>
if ~exist('runnums', 'var') || isempty(runnums)
    
    runnums = cell(1,length(tasks));
    for ii = 1:length(tasks)
        files = dir(fullfile(sessionDir, 'func', sprintf('*task-%s_*bold.nii*', tasks{ii})));
        modality = 'mri';
        if isempty (files)
            files = dir(fullfile(sessionDir, 'meg', sprintf('*task-%s_*meg.sqd', tasks{ii})));
            modality = 'meg';
            if isempty (files)
                files = dir(fullfile(sessionDir, 'ieeg', sprintf('*task-%s_*ieeg.eeg', tasks{ii})));
                modality = 'ecog';
            end
        end
        if ~isempty(files)
            for jj = 1:length(files)
                %runnum{ii}(jj) = str2double(bidsGet(files(jj).name, 'run'));
                runnums{ii}{jj} = bidsGet(files(jj).name, 'run');
            end
        end
    end
end
if ~iscell(runnums), runnums = {{runnums}}; end
if ~iscell(runnums(1)), runnums = {runnums}; end
end