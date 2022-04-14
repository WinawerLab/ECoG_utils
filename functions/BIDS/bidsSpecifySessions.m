function [sessions] = bidsSpecifySessions(projectDir, subject, sessions)

% Specify tasks and run numbers and verify paths for BIDS session
% [sessions] = bidsSpecifyData(projectDir, subject, [session])
%
% Input
%     projectDir:       path where the BIDS projects lies (string)
%     subject:          BIDS subject name (string, all lower case)
%     sessions:         BIDS session name (string, all lower case)
%                           default: all sessions in subject
%
% Output
%     sessions:          String. See input
%
% Example 1
%     projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual'; 
%     subject           = 'p03';
%     [sessions, tasks, runnum] = bidsSpecifySessions(projectDir, subject)

% <projectDir>
if ~exist('projectDir', 'var') || isempty(projectDir)
    error('projectDir not defined'); 
end    

% <subject>
if ~exist('subject', 'var') || isempty(subject)
    error('subject not defined'); 
end

subjectDir = fullfile(projectDir, sprintf('sub-%s', subject));
if ~exist(subjectDir, 'dir')
    error('subject dir not found: %s. Is the server mapped?', subjectDir); 
end

%% Set the optional inputs

% <session>
if ~exist('sessions', 'var') || isempty(sessions)
    d = dir(fullfile(subjectDir, 'ses-*'));
    if isempty(d), error('No session folders found in %s', subjectDir); end
    
    sessions = cell(1,length(d));
    for ii = 1:length(d)
        sessions{ii} = bidsGet(d(ii).name, 'ses');
    end
end

