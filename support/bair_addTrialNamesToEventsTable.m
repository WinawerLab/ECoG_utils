function events = bair_addTrialNamesToEventsTable(events)

fprintf('[%s] Adding trial names...\n',mfilename);

% update events table 
trial_name = cell(height(events),1);

spat_patterns = {'CRF-1','CRF-2','CRF-3','CRF-4','CRF-5',...
                 'GRATING','PLAID','CIRCULAR',...
                 'SPARSITY-1','SPARSITY-2', 'SPARSITY-3', 'SPARSITY-4'};

temp_patterns = {'ONEPULSE-1','ONEPULSE-2','ONEPULSE-3','ONEPULSE-4','ONEPULSE-5','ONEPULSE-6',...
                 'TWOPULSE-1','TWOPULSE-2','TWOPULSE-3','TWOPULSE-4','TWOPULSE-5','TWOPULSE-6'};

spat_objects = {'FACES','LETTERS' 'HOUSES'};

dot_task = {'TASK', 'REST'};

for ii = 1:height(events)
    
    if events.trial_type(ii) == 1
        trial_name{ii} = 'HRFPATTERN';
    elseif events.trial_type(ii) == 2
        trial_name{ii} = 'CONTRASTREVERSEDHRFPATTERN';
    elseif max(events.trial_type(ii) == 3:106)
        trial_name{ii} = 'PRF';
    elseif max(events.trial_type(ii) == 107:118)
        trial_name{ii} = spat_patterns{events.trial_type(ii) == 107:118};
	elseif max(events.trial_type(ii) == 119:130)
        trial_name{ii} = temp_patterns{events.trial_type(ii) == 119:130};
	elseif max(events.trial_type(ii) == 131:133)
        trial_name{ii} = spat_objects{events.trial_type(ii) == 131:133};
    elseif max(events.trial_type(ii) == [134 135])
        trial_name{ii} = dot_task{events.trial_type(ii) == [134 135]};
    elseif events.trial_type(ii) == 255
        trial_name{ii} = 'BLANK';
	elseif events.trial_type(ii) == 256
        trial_name{ii} = 'TASKONSET';
    end 
end

events.trial_name = trial_name;

end