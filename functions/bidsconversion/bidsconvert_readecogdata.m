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
else
    session_inx = 1;
end
fileName = [dataFiles(session_inx).folder filesep dataFiles(session_inx).name];    
fprintf('[%s] Reading %s ...\n', mfilename, fileName);
hdr = ft_read_header(fileName);
data = ft_read_data(fileName);
% To read in EDF data with channels with different sampling rates: data = edf2fieldtrip(fileName);
end