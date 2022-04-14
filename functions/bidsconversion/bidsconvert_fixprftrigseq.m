function [trigSeq] = bidsconvert_fixprftrigseq(stimulus)

requestedtriggers = find(stimulus.trigSeq);
expected = round(linspace(requestedtriggers(1), requestedtriggers(end-1), 224));

% sample_rate = 512;
% for ii = 1:224
%     [gap, idx] = min(abs(expected(ii)-requestedtriggeronsets));
% 
%     if gap < .100 * sample_rate
%         corrected(ii) = requestedtriggeronsets(idx);
%     else
%         corrected(ii) = expected(ii);
%     end
% end
% corrected = expected;

% figure; hold on
% stem(requestedtriggers, ones(1, length(requestedtriggers))); hold on;
% stem(expected, ones(1,224)*1.1);
fprintf('[%s] Creating new PRF triggers by linear interpolation \n', mfilename)


trigSeq = zeros(size(stimulus.trigSeq));
%trigSeq = corrected;
trigSeq(expected) = 1;
        
end