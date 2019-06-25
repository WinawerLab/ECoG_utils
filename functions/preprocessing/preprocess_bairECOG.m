function preprocess_bairECOG(dataPth, sub_label, ses_label, specs)

% This function preprocesses the BIDS-formatted BAIR ECoG data. 

% This is the order of operations:
%
% [1] Reading in data, events and channel info
% [2] Computing matches with visual atlases (wang and benson)
% [3] Common average reference (regression based - uses 'good' channels)
% [4] Broadband computation
% [5] Segmentation (epoching) (includes definition of blank 'trials')
% [6] Baseline correction (optional)
% [7] Computation of spectra
% 
% OUTPUT:
% - data: struct with entire raw and broadband time courses. 
% - trials: struct with evoked, broadband, spectra for each trial. 

%% [0] Path specifications

% Output paths specs
dataDir = fullfile(dataPth, sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label), 'ieeg');
saveDir = fullfile(dataPth, 'derivatives', 'preprocessed', sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label));

% Check whether we have write directories
if ~exist(saveDir, 'dir'); mkdir(saveDir);end

%% [1] Read in ECoG data

fprintf('[%s] Reading data for sub-%s, ses-%s...\n',mfilename, sub_label, ses_label);
[ftdata, events, channels]  = ecog_readBIDSData(dataDir, sub_label, ses_label);

% Check if there are trial_names, if not, add them
if  max(contains(events.Properties.VariableNames, 'trial_name')) == 0
    fprintf('[%s] Events file lacks trial names - adding them now .\n',mfilename);
    events = bair_addTrialNamesToEventsTable(events);
end

% Resample and shift the UMCU data 
if contains(sub_label, 'umcu')
    
    % Resample
    cfg = [];
    cfg.resamplefs      = 512; %frequency at which the data will be resampled (default = 256 Hz)
    cfg.detrend         = 'no';%'no' or 'yes', detrend the data prior to resampling (no default specified, see below)
    %cfg.feedback        = 'no', 'text', 'textbar', 'gui' (default = 'text')
    [ftdata] = ft_resampledata(cfg, ftdata);
    ftdata.hdr.Fs = ftdata.fsample;
    
    % Shift 
    shiftInSeconds = 0.062; % 62 ms
	shiftInSamples = round(shiftInSeconds/(1/512)); % assuming sample rate of 512
    events.onset = events.onset + shiftInSeconds;
    %events.event_sample = events.event_sample + shiftInSamples; 
end

%% [1] Remove non-relevant channels

% % Remove redundant channels (e.g. DC channels, channels without coordinates)
% relevant_channels = find(contains(lower(channels.type), {'ecog', 'seeg','ecg','trig'}));
% 
% ftdata.trial{1} = ftdata.trial{1}(relevant_channels,:);
% ftdata.label = ftdata.label(relevant_channels);
% ftdata.hdr.label = ftdata.hdr.label(relevant_channels);
% ftdata.hdr.chantype = ftdata.hdr.chantype(relevant_channels);
% ftdata.hdr.chanunit = ftdata.hdr.chanunit(relevant_channels);
% 
% channels = channels(relevant_channels,:);

%% [2] Compute matches with visual regions

E2NSpecs = [];
E2NSpecs.pID           = sub_label; % patient ID number
E2NSpecs.patientPool   = 'BAIR';
switch specs.make_plots
    case 'yes'
        E2NSpecs.plotmesh      = 'both';
        E2NSpecs.plotelecs     = 'yes';
    case 'no'
        E2NSpecs.plotmesh      = 'none';
        E2NSpecs.plotelecs     = 'no';
end
visualelectrodes       = electrode_to_nearest_node(E2NSpecs, dataDir);

%% [3] Do a common average reference using the good channels

% NYU preprocessed BIDS data contain a status column, but umcu BIDS data 
% does not; in that case, use all channels as reference... 
if max(contains(channels.Properties.VariableNames, 'status'))>0
    good_channels = find(contains(channels.status, 'good'));
else
    good_channels = find(contains(lower(channels.type),{'ecog','seeg'}));
end
signal = ecog_carRegress(ftdata.trial{1}, good_channels);

%% [3] DIAGNOSTICS: Look at the effect of CAR

switch specs.make_plots
    case 'yes'
        
        figure('Name', 'car regress');
        subplot(1,2,1); hold on
        channel_plot = good_channels(1);
        plot(ftdata.time{1},ftdata.trial{1}(channel_plot,:),'k')
        plot(ftdata.time{1},signal(channel_plot,:),'g')
        legend({'before CAR','after CAR'})
        xlabel('Time (s)'); ylabel('Voltage');
        title(ftdata.label(channel_plot));

        subplot(1,2,2),hold on
        [pxx,freqs] = pwelch(ftdata.trial{1}(good_channels,:)',ftdata.fsample,0,ftdata.fsample,ftdata.fsample);
        [pxx2,~] = pwelch(signal',ftdata.fsample,0,ftdata.fsample,ftdata.fsample);
        plot(freqs,pxx(:,channel_plot),'k')
        plot(freqs,pxx2(:,channel_plot),'g'); 
        set(gca, 'YScale', 'log')
        xlabel('Frequency (Hz)'); ylabel('Log power');
        title(ftdata.label(channel_plot));
        
        % save
end

%% [4] Compute time-varying broadband 

[broadband, methodstr] = ecog_extractBroadband(signal', ftdata.fsample, specs.bb_method, specs.bb_bands);

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
        
        % check for a trigger channel
        if max(contains(data.channels.type,'trig'))
            eltomatch = data.channels.name{find(contains(data.channels.type,'trig'))};
        elseif max(contains(data.channels.name, 'MKR')) % Utrecht data
            eltomatch = data.channels.name{find(contains(data.channels.name,'MKR'))};
        else
            fprintf('[%s] Could not find trigger channel, not plotting triggers\n', mfilename)
            eltomatch = [];
        end
        if ~isempty(eltomatch)
            el = ecog_matchChannels(eltomatch, data);
            ecog_plotFullTimeCourse(data,'raw', el); title('triggers');
            set(gcf, 'Name', 'triggers');
        end
        
        % save plot
end

% % Pick one or more of the electrodes with coverage in visual areas: 
% disp(data.viselec.benson14_varea);

switch specs.make_plots
     case 'yes'

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
blank_trials.events.trial_name = repmat('BLANK', [height(blank_trials.events) 1]);
blank_trials.events.onset = blank_trials.events.onset-blank_epoch(1);
blank_trials.events.duration = repmat(blank_epoch(2)-blank_epoch(1),[height(blank_trials.events) 1]); 
if isfield(summary(trials.events), 'ISI'), blank_trials.events.ISI=[]; end
if isfield(summary(trials.events), 'event_sample')
    blank_trials.events.event_sample = blank_trials.events.event_sample - round(blank_epoch(2)-blank_epoch(1)*data.hdr.Fs);
end

% %% [6] Compute spectra per trial
% 
% epochs = zeros(size(trials.evoked,1), size(trials.evoked,3), size(trials.evoked,2));
% 
% for ii = 1:size(trials.evoked,3)
%     epochs(:,ii,:) = trials.evoked(:,:,ii);
% end
% 
% % dora settings
% fft_w = window(@hann,200); % window width
% fft_ov = 100; % overlap
% % do not regress ERP here, because regressing out average evoked response with only a few trials can hurt 
% reg_erp = 0; 
% fft_t = trials.time>0 & trials.time<=.5; % time segment for spectrum
% 
% [f,spectra] = ecog_spectra(epochs,stims,fft_w,fft_t,fft_ov,trials.fsample,reg_erp);

%% Save preprocessed and epoched data for further analysis
 
% fprintf('[%s] Saving preprocessed continuous data...\n',mfilename);
% %saveName = fullfile(saveDir, sprintf('sub-%s_ses-%s_preproc_bbmethod%d_bandwidth%d', sub_label, ses_label, bbmethod, bands(1,2)-bands(1,1)));
% saveName = fullfile(saveDir, sprintf('sub-%s_ses-%s_preproc.mat', sub_label, ses_label));
% save(saveName, 'data'); 

%saveName = fullfile(saveDir, sprintf('sub-%s_ses-%s_epoched_bbmethod%d_bandwidth%d', sub_label, ses_label, bbmethod, bands(1,2)-bands(1,1)));
saveName = fullfile(saveDir, sprintf('sub-%s_ses-%s_epoched.mat', sub_label, ses_label));

fprintf('[%s] Saving preprocessed epoched data...\n',mfilename);
if exist(saveName,'file')
    fprintf('[%s] Warning: epoched data file already exists!\n',mfilename);
    switch specs.overwrite
        case 'yes'
            fprintf('[%s] Overwriting...\n',mfilename);
        case 'no'
            fprintf('[%s] Saving as separate file with timestamp...\n',mfilename);
            saveName = fullfile(saveDir, sprintf('sub-%s_ses-%s_epoched_%s', sub_label, ses_label, datestr(now,30)));
    end
end

save(saveName, 'trials', 'blank_trials', '-v7.3'); 
fprintf('[%s] Finished preprocessing!\n',mfilename);


