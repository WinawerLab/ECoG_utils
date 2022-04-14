function [outliers,pxx,f] = ecog_plotChannelSpectra(data, chan_select, hdr)

    
if nargin < 3 || isempty(hdr)
    if isstruct(data)
        % this means data was read in with ft_preprocessing instead of
        % ft_read_data
        hdr = data.hdr;
    else
        error('need header info');
    end
end

if isstruct(data)
    % this means data was read in with ft_preprocessing instead of
    % ft_read_data
    data = data.trial{1};    
end

% This is one way to calculate the power spectrum and make a plot, the nice
% thing is that when you click on the lines it returns the channel numbers
figure;
spectopo(data(chan_select,:),size(data,2),hdr.Fs);

% This is another way to look at the power spectrum and make a plot. This
% returns the frequency (f) and power (pxx) and you can check for outliers.
[pxx,f] = pwelch(data(chan_select,:)',hdr.Fs,0,hdr.Fs,hdr.Fs);
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
    disp(chan_select(idx));
    disp(hdr.label(chan_select(idx)));
    outliers = [outliers chan_select(idx)];
end

idx = find(mn<(mean(mn,2)-2*sd));
if ~isempty(idx)
    disp('Two SD below the mean:');
    disp(chan_select(idx));
    disp(hdr.label(chan_select(idx)));
    outliers = [outliers chan_select(idx)];
end

% idx = find(mn>(mean(mn,2)+2*sd));
% if ~isempty(idx)
%     disp('Two SD above the mean:');
%     disp(idx);
%     disp(hdr.label(idx));
%     outliers = [outliers idx];
% end
% 
% idx = find(mn<(mean(mn,2)-2*sd));
% if ~isempty(idx)
%     disp('Two SD below the mean:');
%     disp(idx);
%     disp(hdr.label(idx));
%     outliers = [outliers idx];
% end

end

