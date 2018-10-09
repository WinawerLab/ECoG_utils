function [outliers,pxx,f] = ecog_plotChannelSpectra(data, chan_select)


% This is one way to calculate the power spectrum and make a plot, the nice
% thing is that when you click on the lines it returns the channel numbers
figure;
spectopo(data.trial{1}(chan_select,:),size(data.trial{1},2),data.fsample);

% This is another way to look at the power spectum and make a plot. This
% returns the frequency (f) and power (pxx) and you can check for outliers.
[pxx,f] = pwelch(data.trial{1}(chan_select,:)',data.fsample,0,data.fsample,data.fsample);
%figure,plot(f,log10(pxx));
%xlabel('Frequency (Hz)'); ylabel('Log power');

% Identify channels whose average logpower across frequencies is 2 sd
% above/below the mean across channels:
mn = mean(log10(pxx),1);
sd = std(mn,0,2);

outliers = [];

idx = find(mn>(mean(mn,2)+2*sd));
if ~isempty(idx)
    disp('Two SD above the mean:');
    disp(idx);
    disp(data.label(idx));
    outliers = [outliers idx];
end

idx = find(mn<(mean(mn,2)-2*sd));
if ~isempty(idx)
    disp('Two SD below the mean:');
    disp(idx);
    disp(data.label(idx));
    outliers = [outliers idx];
end

end

