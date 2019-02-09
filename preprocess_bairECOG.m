function preprocess_bairECOG(dataPth, sub_label, ses_label, specs)

% This function preprocesses the BIDS-formatted BAIR ECoG data. 

% This is the order of operations:
%
% [1] Reading in data, events and channel info
% [2] Computing matches with visual atlases (wang and benson)
% [3] Common average reference (regression based)
% [4] Broadband computation
% [5] Segmentation (epoching)
% 
% OUTPUT:
% - data: struct with entire raw and broadband time courses. 
% - trials: struct with evoked and a broadband response for each stimulus. 

%% [0] Path specifications

% Output paths specs
dataDir = fullfile(dataPth, sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label), 'ieeg');
saveDir = fullfile(dataPth, 'derivatives', 'preprocessed', sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label));

% Check whether we have write directories
if ~exist(saveDir, 'dir'); mkdir(saveDir);end

%% [1] Read in ECoG data

fprintf('[%s] Reading data...\n',mfilename);
[ftdata, events, channels]  = ecog_readBIDSData(dataDir, sub_label, ses_label);

% Check if there are trial_names, if not, add them
%if ~isfield(summary(events), 'trial_name')
if  max(contains(events.Properties.VariableNames, 'trial_name'))>0
    fprintf('[%s] Events file lacks trial names - adding them now .\n',mfilename);
    events = bair_addTrialNamesToEventsTable(events);
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

% NYU preprocessed data contain a status column, but umcu data does not;
% in that case, use all channels as reference... (check w Gio whether he
% removed bad channels)
%if isfield(summary(channels), 'status')
if max(contains(channels.Properties.VariableNames, 'status'))>0
    good_channels = find(contains(channels.status, 'good'));
else
    good_channels = find(contains(lower(channels.type),{'ecog','seeg'}));
end
signal = ecog_carRegress_dora(ftdata.trial{1}, good_channels);

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
        [pxx,freqs] = pwelch(ftdata.trial{1}(contains(lower(channels.type),{'ecog','seeg'}),:)',ftdata.fsample,0,ftdata.fsample,ftdata.fsample);
        [pxx2,~] = pwelch(signal',ftdata.fsample,0,ftdata.fsample,ftdata.fsample);
        plot(freqs,pxx(:,channel_plot),'k')
        plot(freqs,pxx2(:,channel_plot),'g'); 
        set(gca, 'YScale', 'log')
        xlabel('Frequency (Hz)'); ylabel('Log power');
        title(ftdata.label(channel_plot));
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
end

% % Pick one or more of the electrodes with coverage in visual areas: 
% disp(data.viselec.benson14_varea);

switch specs.make_plots
     case 'yes'

        if exist('visualelectrodes', 'var') % % e.g.
            eltomatch = visualelectrodes.benson14_varea.elec_labels(1); %MO01-04
        else
            eltomatch = data.channels.name(1);
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

%% Save preprocessed and epoched data for further analysis
% 
% disp('saving preprocessed data...');
% saveName = fullfile(saveDir, sprintf('sub-%s_ses-%s_preproc_bbmethod%d_bandwidth%d', sub_label, ses_label, bbmethod, bands(1,2)-bands(1,1)));
% save(saveName, 'data'); 

%saveName = fullfile(saveDir, sprintf('sub-%s_ses-%s_epoched_bbmethod%d_bandwidth%d', sub_label, ses_label, bbmethod, bands(1,2)-bands(1,1)));
saveName = fullfile(saveDir, sprintf('sub-%s_ses-%s_epoched.mat', sub_label, ses_label));

fprintf('[%s] Saving epoched data...\n',mfilename);
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

save(saveName, 'trials', '-v7.3'); 
fprintf('[%s] Finished preprocessing!\n',mfilename);


