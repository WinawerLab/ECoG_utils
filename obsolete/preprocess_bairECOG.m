function preprocess_bairECOG(dataPth, sub_label, ses_label, specs)

% This function preprocesses the BIDS-formatted BAIR ECoG data. 

% This is the order of operations:
%
% [1] Reading in data, events and channel info
% [2] Computing matches with visual atlases (wang and benson)
% [3] Common average reference (regression based - uses 'good' channels);
%     performed separately for depth, grids, and other surface elecs, 
%     and if there are multiple grids, separately for separate grids
% [4] Broadband computation
% [5] Segmentation (epoching) (includes definition of blank 'trials')
% [6] Save out epoched data to 'derivatives/preprocessed' BIDS folder
% 
% OUTPUT:
% - data: struct with entire raw and broadband time courses. 
% - trials: struct with evoked, broadband, spectra for each trial. 

%% [0] Path specifications

% Output paths specs
dataDir = fullfile(dataPth, sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label), 'ieeg');
saveDir = fullfile(dataPth, 'derivatives', 'preprocessed', sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label));
figSaveDir = fullfile(dataPth, 'derivatives', 'preprocessed', sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label), 'figures', 'preprocessing');

% Check whether we have write directories
if ~exist(saveDir, 'dir'); mkdir(saveDir);end
if ~exist(figSaveDir, 'dir'); mkdir(figSaveDir);end

%% [1] Read in ECoG data

fprintf('[%s] Reading data for sub-%s, ses-%s...\n',mfilename, sub_label, ses_label);
[ftdata, events, channels]  = ecog_readBIDSData(dataDir, sub_label, ses_label);

fprintf('[%s] Checking for trial names...\n',mfilename);

% Check if there are trial_names, if not, add them
if  max(contains(events.Properties.VariableNames, 'trial_name')) == 0
    fprintf('[%s] Events file lacks trial names - adding them now .\n',mfilename);
    events = bair_addTrialNamesToEventsTable(events);
end

% Resample and SHIFT the UMCU data 
if contains(sub_label, {'chaam', 'intraop'})
    fprintf('[%s] This is a umcu patient: resampling the data\n',mfilename);
    % Resample
    cfg = [];
    cfg.resamplefs      = 512;  % frequency at which the data will be resampled 
    cfg.detrend         = 'no'; %'no' or 'yes', detrend the data prior to resampling (no default specified, see below)
    [ftdata] = ft_resampledata(cfg, ftdata);
    ftdata.hdr.Fs = ftdata.fsample;
end
if contains(sub_label, 'chaam')
    fprintf('[%s] This is umcu patient chaam: shifting the data\n',mfilename);
    % Shift 
    shiftInSeconds = 0.062; % 62 ms
	shiftInSamples = round(shiftInSeconds/(1/512)); % assuming sample rate of 512
    events.onset = events.onset + shiftInSeconds;
    if isfield(summary(events), 'event_sample')
        events.event_sample = events.event_sample + shiftInSamples; 
    end
end
if contains(sub_label, 'intraop')
    channels.status = repmat({'good'}, [height(channels) 1]);
    channels.status(80:96,:) = {'bad'};
end
%% [2] Compute matches with visual regions
fprintf('[%s] Computing matches with visual atlases...\n',mfilename);
E2NSpecs = [];
E2NSpecs.pID           = sub_label; % patient ID number
switch specs.make_plots
    case 'yes'
        close all;
        E2NSpecs.plotmesh      = 'both';
        E2NSpecs.plotelecs     = 'yes';
    case 'no'
        E2NSpecs.plotmesh      = 'none';
        E2NSpecs.plotelecs     = 'no';
end
visualelectrodes       = electrode_to_nearest_node(E2NSpecs, dataDir);

% Add visual area names (W and B) ecc, angle, sigma to channels table
[channels] = bair_addVisualAtlasNamesToChannelTable(channels,visualelectrodes);

switch specs.make_plots
    case 'yes'
        fprintf('[%s] Saving pial mesh figures...\n',mfilename);
        
        fig = get(groot,'CurrentFigure');
        if ~isempty(fig)
            nPlots = get(gcf,'Number');
            for ii = 1:nPlots
                atlasName = get(ii,'Name');
                atlasName = strsplit(atlasName);
                if length(atlasName) == 1
                    atlasName = 'none';
                else
                    atlasName = atlasName{2};
                end
                plotName = sprintf('pialmesh_atlas-%s',atlasName);
                saveas(ii, fullfile(figSaveDir, sprintf('%s-%s-%s',sub_label, ses_label,plotName)), 'epsc');
            end
        end
end
close all;
%% [3] Do a common average reference using the good channels

fprintf('[%s] Performing Common Average reference...\n',mfilename);

[signal, channels, INX, INXNames] = ecog_performCAR(ftdata.trial{1}, channels);

% DIAGNOSTICS: Look at the effect of CAR
switch specs.make_plots
    case 'yes'

    fprintf('[%s] Saving CAR figures...\n',mfilename);

    for ii = 1:length(INX)
    
        chan_index = INX{ii};
    
        if ~isempty(chan_index)
        
            good_channels = find(contains(channels(chan_index,:).status, 'good'));
            channel_plot = chan_index(good_channels(1));

            figure('Name', sprintf('car regress %s', INXNames{ii}));
            
            % Plot the time courses before and after CAR
            subplot(1,2,1); hold on
            plot(ftdata.time{1},ftdata.trial{1}(channel_plot,:),'k')
            plot(ftdata.time{1},signal(channel_plot,:),'g')
            legend({'before CAR','after CAR'})
            xlabel('Time (s)'); ylabel('Voltage');set(gca, 'FontSize', 16);
            title(ftdata.label(channel_plot));

            % Plot the spectra before and after CAR
            subplot(1,2,2),hold on
            [pxx,freqs] = pwelch(ftdata.trial{1}(channel_plot,:)',ftdata.fsample,0,ftdata.fsample,ftdata.fsample);
            [pxx2,~] = pwelch(signal(channel_plot,:)',ftdata.fsample,0,ftdata.fsample,ftdata.fsample);
            plot(freqs,pxx,'k')
            plot(freqs,pxx2,'g'); 
            set(gca, 'YScale', 'log')
            xlabel('Frequency (Hz)'); ylabel('Log power');set(gca, 'FontSize', 16);
            title(ftdata.label(channel_plot));
            
            set(gcf, 'Position',[1000 700 1100 600]); 
            saveas(gcf, fullfile(figSaveDir, sprintf('%s-%s-CAR-%s',sub_label, ses_label, INXNames{ii})), 'epsc');close;
        end
    end
end

%% [4] Compute time-varying broadband 

[broadband, methodstr] = ecog_extractBroadband(signal', ftdata.fsample, specs.bb_method, specs.bb_bands);
broadband = broadband';

% Create a data structure to save
data = struct();

data.hdr        = ftdata.hdr;
data.time       = ftdata.time{1};
data.raw        = ftdata.trial{1}; 
data.car_reref  = signal;
data.broadband  = broadband;
data.bb_bands   = specs.bb_bands;
data.bb_method  = methodstr;
data.events     = events;
data.channels   = channels;
data.viselec    = visualelectrodes;
data.cfg        = ftdata.cfg;

%% [4] DIAGNOSTICS: Plot filtered time course and broadband for channel(s) of interest

% Check whether the event onsets align with the triggers by plotting
% trigger channel itself Note that there may be slight offsets because we
% segment on the flips, not the triggers!

switch specs.make_plots
    case 'yes'
        
        fprintf('[%s] Plotting the trigger channel...\n',mfilename);
        % check for a trigger channel
        if max(contains(data.channels.type,'trig'))
            eltomatch = data.channels.name{find(contains(data.channels.type,'trig'))};
        elseif max(contains(data.channels.name, 'MKR')) % Utrecht data
            eltomatch = data.channels.name{find(contains(data.channels.name,'MKR'))};
        else
            fprintf('[%s] Could not find trigger channel, not plotting triggers\n', mfilename);
            eltomatch = [];
        end
        if ~isempty(eltomatch)
            el = ecog_matchChannels(eltomatch, data);
            ecog_plotFullTimeCourse(data,'raw', el); title('triggers'); 
            set(gcf, 'Name', 'triggers');
        end
        
        % save plot
        saveas(gcf, fullfile(figSaveDir, sprintf('%s-%s-triggerchannel',sub_label, ses_label)), 'epsc'); close

end

% % Pick one or more of the electrodes with coverage in visual areas: 
% disp(data.viselec.benson14_varea);

switch specs.make_plots
     case 'yes'
        fprintf('[%s] Plotting a data channel...\n',mfilename);

        eltomatch = data.channels.name(1);
        if exist('visualelectrodes', 'var') && ~ isempty(visualelectrodes)% % e.g.
            if ~isempty(visualelectrodes.benson14_varea)
                if ~isempty(visualelectrodes.benson14_varea.elec_labels)
                    eltomatch = visualelectrodes.benson14_varea.elec_labels(1); %MO01-04
                end
            end
        end
         
        % % Find matching channel number
        el = ecog_matchChannels(eltomatch, data);
         
        % % Plot voltage and broadband
        % ecog_plotFullTimeCourse(data,'car_reref', el);
        ecog_plotFullTimeCourse(data,'broadband', el); % no smoothing
        % ecog_plotFullTimeCourse(data,'broadband', el, data.hdr.Fs/2); % with smoothing
        set(gcf, 'Name', 'example broadband time course');
        saveas(gcf, fullfile(figSaveDir, sprintf('%s-%s-examplechannel',sub_label, ses_label)), 'epsc'); close;

end

%% [5] Segment the data (all channels, all tasks)

% Epoch the stimulus trials
trials = ecog_epochData(data, specs.epoch);

% Add 'no stimulation' blank trials for non-prf tasks
fprintf('[%s] Constructing blank events...\n',mfilename);

mockdata = data;
if isfield(summary(mockdata.events), 'task_name')
    mockdata.events = mockdata.events(contains(data.events.task_name, {'hrfpattern','spatialpattern', 'temporalpattern', 'spatialobject', 'soc'}),:);
else
    mockdata.events = mockdata.events(~contains(data.events.trial_name, {'PRF', 'BLANK'}),:);
end

% Define an epoch in which there were no stimuli
% hrf: ITI between 3 and 24 seconds, stimduration 0.2s
% spatpat/spatobj/temppat: ITI between 1.25 and 1.75 seconds, stimduration 0.5s

blank_epoch = [specs.epoch(1)-0.5 specs.epoch(1)];
if blank_epoch(1) < -.75
    fprintf('[%s] Warning: blank epoch overlaps with stimulus epoch for spatialpattern/temporalpattern/spatialobject! \n',mfilename);
end

% Epoch the baseline 'trials'
blank_trials = ecog_epochData(mockdata, blank_epoch);
blank_trials.events.trial_name = repmat({'BLANK'}, [height(blank_trials.events) 1]);
blank_trials.events.onset = blank_trials.events.onset-blank_epoch(1);
blank_trials.events.duration = repmat(blank_epoch(2)-blank_epoch(1),[height(blank_trials.events) 1]); 
if isfield(summary(trials.events), 'ISI'), blank_trials.events.ISI= zeros([height(blank_trials.events) 1]); end
if isfield(summary(trials.events), 'event_sample')
    blank_trials.events.event_sample = blank_trials.events.event_sample - round(blank_epoch(2)-blank_epoch(1)*data.hdr.Fs);
end

%% [6] Save epoched data for further analysis
 
saveName = fullfile(saveDir, sprintf('sub-%s_ses-%s_epoched.mat', sub_label, ses_label));

fprintf('[%s] Saving preprocessed epoched data to %s...\n',mfilename, saveName);
if exist(saveName,'file')
    warning('[%s] Epoched data file already exists!',mfilename);
    switch specs.overwrite
        case 'yes'
            fprintf('[%s] Overwriting...\n',mfilename);
        case 'no'
            saveName = fullfile(saveDir, sprintf('sub-%s_ses-%s_epoched_%s', sub_label, ses_label, datestr(now,30)));
            fprintf('[%s] Saving preprocessed epoched data to %s...\n',mfilename, saveName);
    end
end

save(saveName, 'trials', 'blank_trials', '-v7.3'); 
fprintf('[%s] Finished preprocessing!\n',mfilename);


