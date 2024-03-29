function [avg_data, events, n_avg, events_idx] = ecog_averageEvents(data,events,average,avg_fun)

% Description: 
%
% [data,events,n_avg,events_idx] = ECOG_AVERAGEEVENTS(data,events,average_type,[average_method])
% 
% takes average across events based on average_type.
% 
% Input
% - data           = numeric matrix (t x events x channels)
% - events         = events table
% - average_type   = 'seesions', 'runs', 'stimuli', 'trials', 'none'
%                    'stimuliINsessions', 'stimuliINruns'
%                    'trialsINsessions', 'trialsINruns'
% - average_method = average method (default = @mean)
% 
% Output
% - n_avg          = numbers of average (also outpus in events)
% - events_idx     = cell-array of index of averaged trials

% Dependency: SetDefault

% 20220803 Yuasa

%% parameter setting
narginchk(3,4);
SetDefault('avg_fun',@mean);

%-- check data validity
datsiz = size(data);
ntrl   = datsiz(2);
assert(height(events) == ntrl, 'The numbers of events does not match input data');

%% average
%-- prepare for average
if ~isfield(summary(events),'stim_file_index')
    events_cp = events;
    events_cp.stim_file_index(:) = 0;
    events_cp.stim_orig_idx      = [1:height(events_cp)]';
    events_cp = sortrows(events_cp,{'task_name','session_name','run_name','trial_type'});
    avg_group = findgroups(events(:,{'task_name','session_name','run_name'}));
    for ii=reshape(unique(avg_group),1,[])
        events_cp.stim_file_index(avg_group==ii) = 1:sum(avg_group==ii);
    end
    events_cp = sortrows(events_cp,'stim_orig_idx');
    events.stim_file_index = events_cp.stim_file_index;
end
switch average
    case {'none'},      avg_group = [1:ntrl]';
    case {'sessions'},  avg_group = findgroups(events(:,{'task_name','run_name','stim_file_index'}));
    case {'runs'},      avg_group = findgroups(events(:,{'task_name','stim_file_index'}));
    case {'stimuli'},   avg_group = findgroups(events(:,{'task_name','trial_type'}));
    case {'trials'},    avg_group = findgroups(events.task_name);
    case {'all'},       avg_group = ones(ntrl,1);
    case {'stimuliINsessions'}, avg_group = findgroups(events(:,{'task_name','session_name','trial_type'}));
    case {'stimuliINruns'},     avg_group = findgroups(events(:,{'task_name','session_name','run_name','trial_type'}));
    case {'trialsINsessions'},  avg_group = findgroups(events(:,{'task_name','session_name'}));
    case {'trialsINruns'},      avg_group = findgroups(events(:,{'task_name','session_name','run_name'}));
    otherwise,          error('''%s'' is unknown average type',average);
end
if startsWith(average,'trials')
    %-- segregate BLANK
    bslIndex  = contains(events.trial_name, 'BLANK');
    avg_group(bslIndex) = avg_group(bslIndex) + max(avg_group);
    
    %-- modify events
    events.trial_type(~bslIndex) = min(events.trial_type(~bslIndex));
    events.trial_name(~bslIndex) = events.task_name(~bslIndex);
end
if matches(average,'all')
    %-- modify events
    events.trial_type(:) = min(events.trial_type(:));
    if iscell(events.trial_name)
        events.trial_name(:) = {'ALL'};
    else
        events.trial_name(:) = 'ALL';
    end
end
    
n_avg     = groupcounts(avg_group);
if ~all(diff(n_avg)==0),    warning('[%s] %s not consist with the same time series',mfilename,average);  end

%-- take average across repeats (t x events x chan -> t x averaged events x chan)
datsiz(2)   = length(n_avg);
avg_data    = zeros(datsiz);
events_idx  = cell(length(n_avg),1);
for iavg = 1:length(n_avg)
    avg_data(:,iavg,:,:)   = avg_fun(data(:,avg_group==iavg,:,:),2,'omitnan');
    events_idx{iavg}       = find(avg_group==iavg);
end

events       = events(cellfun(@(I) I(1),events_idx),:);
events.n_avg = n_avg;

end
