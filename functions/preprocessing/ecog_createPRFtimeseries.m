function ts = ecog_createPRFtimeseries(epochs,t, time_win, stimInx, normalizeEpochs)


% Determine number of runs and stimuli for this subject
[~, nTrials, nChans] = size(epochs);
nStim = length(stimInx);    
nRuns = nTrials/nStim; 

% Normalize epochs, within run?
if normalizeEpochs
    run_indices = nan(nStim,nRuns); 
    for jj = 1:nRuns
        run_indices(:,jj) = ones(nStim,1)*jj;
    end
    [epochs] = ecog_normalizeEpochs(epochs, t, [], [], run_indices);
end

% Compute average broadband response in time window
trials = squeeze(mean(epochs(t>time_win(1) & t<time_win(2),:,:),1)); 

% Transpose to have channels in first dimension
trials = trials';

% Reshape to separate individual runs    
ts = reshape(trials,[nChans nStim nRuns]);
    
    
end