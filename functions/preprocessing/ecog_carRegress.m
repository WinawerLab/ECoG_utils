function [signal] = ecog_carRegress(signal, chans2incl)
% common average referencing function with regressing out the mean, rather
% than subtracting out the mean
%
% dh - Oct 2010
% signal = electrodes X samples
% 
% updates IG Aug 2019: CAR is applied to all channels (including bad ones),
% but only the good channels go into the average. 

if size(signal,2) < size(signal,1) % signal samples X electrodes
    disp('transpose signal to be electrodes X samples')
    return
end

ca_signal = mean(signal(chans2incl,:),1);

% regress off the mean signal
for kk = 1:size(signal,1) % elecs
    %disp(['elec ' int2str(kk)])
    %if ismember(kk,chans2incl)
        [B,BINT,R] = regress(signal(kk,:)',ca_signal');
        signal(kk,:)=R';
    %end
end

