function [data, hdr] = bidsconvert_readecogdata(dataReadDir, ses_label)

if nargin < 2 || isempty(ses_label)
    ses_label = 'nyuecog01';
end

% READ IN relevant data files %%%%%%%%%%%%%%%%%%

%   - ECoG data file

% Read ECoG data
dataFiles = dir(fullfile(dataReadDir, '*_512*.EDF'));
if length(dataFiles) > 1 
    warning('[%s] Multiple datafiles found: comparing with required ses_label %s', mfilename, ses_label);
    % Figure out which session
    out = regexp(ses_label, '\d*', 'match');
    session_inx = str2double(out{1}); 
    file_inx = [];
    for ii = 1:length(dataFiles)
        if contains(dataFiles(ii).name, sprintf('%d_512.EDF', session_inx)) || contains(dataFiles(ii).name, sprintf('512_%d.EDF', session_inx))
           file_inx = [file_inx ii];
        end
%         if contains(dataFiles(ii).name, num2str(session_inx))
%             file_inx = [file_inx ii];
%         end
    end
    if length(file_inx) ~= 1
        error('[%s] Could not match session name %s to datafile: please check files', mfilename, ses_label);
    end

else
    file_inx = 1;
end
fileName = [dataFiles(file_inx).folder filesep dataFiles(file_inx).name];    
fprintf('[%s] Reading %s ...\n', mfilename, fileName);
hdr = ft_read_header(fileName);
data = ft_read_data(fileName);
% To read in EDF data with channels with different sampling rates: data = edf2fieldtrip(fileName);
end