function [trials] = ecog_epochData(data, epochDur)

fprintf('[%s] Epoching raw and broadband timecourses...\n',mfilename);

% Do this to make sure methods and specs of preproc get inherited by trials
trials = data;

% Remove fields specific to continuous data
if isfield(trials, 'raw'); trials = rmfield(trials, 'raw'); end
if isfield(trials, 'car_reref'); trials = rmfield(trials, 'car_reref'); end
if isfield(trials, 'broadband'); trials = rmfield(trials, 'broadband'); end
if isfield(trials, 'cfg'); trials = rmfield(trials, 'cfg'); end
if isfield(trials, 'hdr'); trials = rmfield(trials, 'hdr'); end

% Determine how many samples to go back and forth to extract epoch
onset_pre = round(epochDur(1)*data.hdr.Fs);
onset_post = round((epochDur(2)*data.hdr.Fs)-1);

%if isfield(summary(data.events), 'event_sample')
if max(contains(data.events.Properties.VariableNames, 'event_sample'))>0
    onsetInx = data.events.event_sample;
else
    onsetInSeconds = data.events.onset;
    onsetInx = round(onsetInSeconds*data.hdr.Fs);
end

% Overwrite trials
trials.time  = epochDur(1):(1/data.hdr.Fs):epochDur(2)-(1/data.hdr.Fs);

% Extract epochs
for ii = 1:length(onsetInx)  
    % broadband
    trials.broadband(:,:,ii) = data.broadband(:,onsetInx(ii)+onset_pre:onsetInx(ii)+onset_post);
    % evoked
    trials.evoked(:,:,ii) = data.car_reref(:,onsetInx(ii)+onset_pre:onsetInx(ii)+onset_post);
    % broadband
    %trials.broadband{ii} = data.broadband(:,onsetInx(ii)+onset_pre:onsetInx(ii)+onset_post);
    % evoked
    %trials.evoked{ii} = data.car_reref(:,onsetInx(ii)+onset_pre:onsetInx(ii)+onset_post);
end

%disp('done');

trials.events   = data.events;
trials.fsample  = data.hdr.Fs;

end