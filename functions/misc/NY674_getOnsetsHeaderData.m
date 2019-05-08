function [onsets, onsets_run, data, hdr, triggerChannels] = NY674_getOnsetsHeaderData()

patientID = 674;
pth = '/Volumes/server/Projects/BAIR/Data/Raw/ECoG/';
datapth = [pth num2str(patientID)];
datafiles = {'NY674_BAIRvisual_Clin1Part1', 'NY674_BAIRvisual_Clin1Part2','NY674_BAIRvisual_Clin2'};

data = [];
for ii = 1:length(datafiles)
    [hdr{ii}, data{ii}] = edfread([datapth filesep datafiles{ii} '.edf']);
end

sample_rate = hdr{1}.frequency(1);
% concatenate Clin1 data in time
data1 = [data{1} data{2}]; 
data2 = data{3};

triggerChannels = [138 68]; % DC10


for tt = 1:2
    switch tt
        case 1
            triggers = data1(triggerChannels(tt),:);  
            triggerind  = [174300 437100 565700 660000 705000 740000 775000 810000 850000 880000 920000 960000 990000 1030000 1065000 1140000 1250000 inf];
        case 2
            triggers = data2(triggerChannels(tt),:);
            triggerind  = [200000 450000 650000 760000 800000 840000 870000 910000 945000 980000 1020000 1055000 1090000 1125000 1160000 1220000 1345000 inf];

    end
    
    triggers = triggers / max(triggers);    
    t = (0:length(triggers)-1)/sample_rate;
    
    
    figure(1); clf
    plot(triggers);
    
    % approximate trigger indices for start of each of 17 experiments, discarding first one
    
    
    
    [~, trigger_onsets] = findpeaks(triggers , 'MinPeakHeight', 0.8, 'MinPeakDistance', .05*sample_rate);
    trigger_onsets = trigger_onsets(trigger_onsets>triggerind(1)); % 2333
    d = [inf diff(trigger_onsets)];
    keep_idx = d > .3*sample_rate;
    trigger_onsets = trigger_onsets(keep_idx);
    
    figure(2), clf
    plot(t, triggers); hold on
    plot(t(trigger_onsets), triggers(trigger_onsets), 'or')
    
    for ii = 1:17
        idx = trigger_onsets>triggerind(ii) & trigger_onsets<=triggerind(ii+1);
        s(ii)=sum(idx);
        trigger_onsets_run{ii} = trigger_onsets(idx);
    end
    
    figure(3), clf
    plot(t, triggers); hold on
    for ii = 1:17
        plot(t(trigger_onsets_run{ii}), triggers(trigger_onsets_run{ii}), 'o')
    end
    
    % fix pRF triggers
    prfRuns = [2 3 16 17];
    
    for r = 1:4
        this_run = prfRuns(r);
        these_onsets = trigger_onsets_run{this_run};
        expected = round(linspace(these_onsets(1), these_onsets(end), 224));
        
        for ii = 1:224
            [gap, idx] = min(abs(expected(ii)-these_onsets));
            
            if gap < .100 * sample_rate
                corrected(ii) = these_onsets(idx);
            else
                corrected(ii) = expected(ii);
            end
        end
        corrected = expected;
        title(r);
        
        figure(4); clf
        stem(these_onsets, ones(1, length(these_onsets))); hold on;
        stem(expected, ones(1,224)*1.1);
        
        
        %pause(1)
        trigger_onsets_run{this_run} = corrected;
    end
    
    
     %% replot
    figure(5), clf
    plot(t, triggers); hold on
    for ii = 1:17
        plot(t(trigger_onsets_run{ii}), ones(size(trigger_onsets_run{ii})), 'o')
    end
    %%
    
    trigger_onsets = [];
    for ii = 1:17
        trigger_onsets = [trigger_onsets trigger_onsets_run{ii}];
    end
    onsets(:,tt) = trigger_onsets;
    onsets_run{tt} = trigger_onsets_run;
end

data{1} = data1;
data{2} = data2;
hdr{2} = hdr{3};
data = data(1:2);
hdr = hdr(1:2);






    