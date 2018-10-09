function [trials] = ecog_epochData(data, epochDur)

trials = struct();

% determine how many samples to go back and forth to extract epoch
onset_pre = epochDur(1)/(1/data.hdr.Fs);
onset_post = epochDur(2)/(1/data.hdr.Fs)-1;

onsetInx = data.events.onset_index;

% Extract epochs
for ii = 1:length(onsetInx)  
    % broadband
    trials.broadband(:,:,ii) = data.broadband(:,onsetInx(ii)+onset_pre:onsetInx(ii)+onset_post);
    % evoked
    trials.evoked(:,:,ii) = data.car_reref(:,onsetInx(ii)+onset_pre:onsetInx(ii)+onset_post);
end

disp('done');

trials.label    = data.hdr.label;
trials.time     = epochDur(1):(1/data.hdr.Fs):epochDur(2)-(1/data.hdr.Fs);
trials.events   = data.events;
trials.fsample  = data.hdr.Fs;

end