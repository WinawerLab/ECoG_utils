function [channels,chanidx] = ecog_HDgridElectrodes(channels,gridHDthresh,visarearate,vischidx)
% [channels,index] = ecog_HDgridElectrodes(channels,[elec_number,vis_coverrate,vis_index])
%   outputs channels which are parts of HDgrids which covers visual area
% 
% Input
%     channels:  	 channel table
%     elec_number:   numbers of electrodes to segrigate HDgrid or not
%                    (default: 64)
%     vis_coverrate: threshold to judge if the HDgrid covers visual area or not
%                    if 0.1, 10% of electrodes need to overlap on visual areas
%                    (default: 0.1)
%     vis_index:     index of channels which are associated to any visual areas
% 
% Output
%     channels:  	channel table
%     index:        logical array to match input and output channels
% 
% See also, ecog_visualElectrodes
% 
% K.Yuasa, 2022

% Atlas fileds
if ~exist('gridHDthresh','var') || isempty(gridHDthresh)
    gridHDthresh = 64;
end
if ~exist('visarearate','var') || isempty(visarearate)
    visarearate  = 0.1;
end
if ~exist('vischidx','var') || isempty(vischidx)
    [~, vischidx] = ecog_visualElectrodes(channels);
end

% Check HDgrid which covered visual area (with over 'visarearate')
chanidx = [];
if istablefield(channels,'group') && any(ismember(channels.group,'HDgrid'))
  chan_grid = contains(channels.group, 'HDgrid');
  gridTypes = string(unique(regexprep(channels.name(chan_grid),'\d*$','')));
  for igrid = reshape(gridTypes,1,[])
      chan_grid = contains(channels.group, 'HDgrid') & startsWith(channels.name,igrid);
      if sum(all([vischidx,chan_grid],2))>(sum(chan_grid)*visarearate)
          chanidx = [chanidx, chan_grid];
      end
  end
else      % if it's unsure channels.group has HDgrid information or not
  chan_grid = startsWith(channels.name, 'GA');
  if sum(chan_grid) > gridHDthresh
      channels.group(chan_grid) = {'HDgrid'};
      if sum(all([vischidx,chan_grid],2))>(sum(chan_grid)*visarearate)
          chanidx = [chanidx, chan_grid];
      end
  else
      channels.group(chan_grid) = {'grid'};
  end
  chan_grid = startsWith(channels.name, 'GB');
  if sum(chan_grid) > gridHDthresh
      channels.group(chan_grid) = {'HDgrid'};
      if sum(all([vischidx,chan_grid],2))>(sum(chan_grid)*visarearate)
          chanidx = [chanidx, chan_grid];
      end
  else
      channels.group(chan_grid) = {'grid'};
  end
  chan_grid = startsWith(channels.name, 'GC');
  if sum(chan_grid) > gridHDthresh
      channels.group(chan_grid) = {'HDgrid'};
      if sum(all([vischidx,chan_grid],2))>(sum(chan_grid)*visarearate)
          chanidx = [chanidx, chan_grid];
      end
  else
      channels.group(chan_grid) = {'grid'};
  end
  chan_grid = startsWith(channels.name, 'G') & ~startsWith(channels.name, {'GA','GB','GC'});
  if sum(chan_grid) > gridHDthresh
      channels.group(chan_grid) = {'HDgrid'};
      if sum(all([vischidx,chan_grid],2))>(sum(chan_grid)*visarearate)
          chanidx = [chanidx, chan_grid];
      end
  else
      channels.group(chan_grid) = {'grid'};
  end
end

% Output
chanidx  = any(chanidx,2);
channels = channels(chanidx,:);

