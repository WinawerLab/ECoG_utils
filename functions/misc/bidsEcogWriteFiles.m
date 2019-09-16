function [fname_out] = bidsEcogWriteFiles(dataPath, subject, session, task, runnum, description, ...
    data, channels, events, ieeg_json, hdr)

% <outputFolder>
if ~exist('description', 'var')
    description = [];
end

fprintf('[%s] Task = %s, Run = %s \n', mfilename, task, runnum);

if ~isempty(description)
    fname_out = sprintf('sub-%s_ses-%s_task-%s_run-%s_desc-%s', subject, session, task, runnum, description);
else
    fname_out = sprintf('sub-%s_ses-%s_task-%s_run-%s', subject, session, task, runnum);
end

sessionDir = fullfile(dataPath, sprintf('sub-%s', subject), sprintf('ses-%s', session));

% Write out the channels file
chanWriteFile = fullfile(sessionDir, 'ieeg', sprintf('%s_channels.tsv', fname_out));
fprintf('[%s] Writing new channels file: %s\n', mfilename, chanWriteFile); 
writetable(channels,chanWriteFile,'FileType','text','Delimiter','\t');
       
% Write out the events file
eventsWriteFile = fullfile(sessionDir, 'ieeg', sprintf('%s_events.tsv', fname_out));
fprintf('[%s] Writing new channels file: %s\n', mfilename, eventsWriteFile); 
writetable(events,eventsWriteFile,'FileType','text','Delimiter','\t');
  
% Save out the data file
dataWriteFile = fullfile(sessionDir, 'ieeg', sprintf('%s_ieeg.eeg', fname_out));
fprintf('[%s] Writing new data file: %s\n', mfilename, dataWriteFile); 
ft_write_data(dataWriteFile, data, 'header', hdr, 'dataformat', 'brainvision_eeg');   

% Save out the ieeg_json file
jsonWriteFile = fullfile(sessionDir, 'ieeg', sprintf('%s_desc-reref_ieeg.json', fname_out));
fprintf('[%s] Writing new ieeg json file: %s\n', mfilename, jsonWriteFile); 
json_options.indent = '    '; 
jsonwrite(jsonWriteFile,ieeg_json, json_options)
                                      
end

