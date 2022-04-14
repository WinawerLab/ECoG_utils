function [] = bidsconvert_anonymizestimulusfiles(stimDir, dry_run)

% Anonymizes BAIR stimulus files according to the subject name provided in
% the file name. Note this assumes the filename has been updated already
% tot the anonymized name!
%
% Will by default run in dry_run mode, i.e. not overwrite any files just
% index the fields. If you want to overwrite the files, input dry_run = 0.
%
% Also updates the year in the params.experimentDateandTime field to 1900,
% to match the data scans files.

% IG BAIR 2022

if ~exist('dry_run', 'var') || isempty(dry_run)
    dry_run = 1;
end

% find all stimulus .mat files in stimDir
fileNames = dir(fullfile(stimDir, '*.mat'));

if dry_run
    fprintf('dry run, not overwriting any files \n');
    answer = {'n'};
else
    prompt = {'y/n'};
    defaults = {'n'};
    answer = inputdlg(prompt, 'You disabled dry run. Are you sure you want to continue?', [1 50], defaults);
end
    
for ii = 1 :length(fileNames)
    
    fileName = fileNames(ii).name;
    subjectName = bidsGet(fileName, 'sub');
    s = load(fullfile(stimDir,fileName));
    fprintf('loading stimFile %s \n', fullfile(stimDir,fileName));
        
    % replace subject names
    if isfield(s, 'fname')
        fprintf('changing fname %s to %s \n', s.fname, fileName);
        s.fname = fileName;
    end
    if isfield(s.params, 'subjID')
        fprintf('changing subjID %s to %s \n', s.params.subjID, subjectName);
        s.params.subjID = subjectName;
    end
    if isfield(s.params, 'pID')
        fprintf('changing pID %s to %s \n', s.params.pID, subjectName);
        s.params.pID = subjectName;
    end
    
    % remove the tempname field, if present
    if isfield(s, 'tempname')
        fprintf('removing tempname %s \n', s.tempname);
        s = rmfield(s, 'tempname');
    end
    if isfield(s, 'D')
        fprintf('removing D field \n');
        s = rmfield(s, 'D');
    end
        
    % update timestamp
    if isfield(s.params, 'experimentDateandTime')
        experimentDateandTimeOld = s.params.experimentDateandTime;
        experimentDateandTimeNew = experimentDateandTimeOld;
        experimentDateandTimeNew(1:4) = '1900';
        fprintf('changing timestamp %s to %s \n', experimentDateandTimeOld, experimentDateandTimeNew);
        s.params.experimentDateandTime = experimentDateandTimeNew;
    end
    
    if strcmpi(answer, 'y')
        save(fullfile(stimDir,fileName), '-struct', 's');
        fprintf('saving stimFile %s \n', fullfile(stimDir,fileName));
    end
    
end

end