function [trials_out,figlist] = ecog_plotGridSpectra(spectra, whichHDgrid, whichTrials, whichTasks, specs)

% [trials_out, p] = ecog_plotGridSpectra(trials, whichHDgrid, whichTrials, whichTasks, specs)
% Plot the entire HD grid
% 
% See also ecog_plotSpectra

% Dependency: ecog_plotSpectra, SetDefault

% 20200306 Yuasa

%% Set options
narginchk(3,inf);
SetDefault('specs.plot.XScale','linear', 0);
SetDefault('specs.plot.YScale','log', 0);
SetDefault('specs.plot.nSubPlots', [], 1);
SetDefault('specs.plot.RotGrid', false, 0);

% whichHDgrid     = upper(whichHDgrid);
%%
%-- correct elecnames (set '%03d')
[eleccat, elecnum] = strtok(spectra.channels.name,int2str(0:9));
for el = 1:length(spectra.channels.name)
    spectra.channels.name{el} = sprintf('%s%03d',eleccat{el},str2double(elecnum{el}));
end

%-- get electrodes list
whichElectrodes = spectra.channels.name(startsWith(spectra.channels.name, whichHDgrid));

%-- set grid parameters
switch whichHDgrid
    case 'GA',  nCol = 8; nRow = 8;
    case 'GB',  nCol = 8; nRow = 16;
end
if specs.plot.RotGrid, tmp = nCol; nCol = nRow; nRow = tmp; end
FullGRID = nCol * nRow;
if isempty(specs.plot.nSubPlots)
    nFig = ceil(nRow / nCol ./ 1.2);
    specs.plot.nSubPlots = [ceil(nRow./nFig), nCol];
end
plRow  = specs.plot.nSubPlots(1);
nFig = ceil(nRow./plRow);

%-- Create a list of electrode names for the overall grid layout
gridList = [];
for ee = 1:FullGRID
    chanName = sprintf('%s%03d',whichHDgrid,ee);
        
    if any(strcmp(chanName,whichElectrodes))
        gridList{ee} = chanName;
    else
        gridList{ee} = [chanName '-nodata'];
    end
end
 
%-- Plot response per electrode, in two separate figures (top and bottom)
v = reshape(1:FullGRID,[nRow,nCol]);
inx = [];
for ee = 1:nFig
    inx{ee} = reshape(v(plRow*(ee-1)+1:min(nRow,plRow*(ee)),:)',1,[]);
end

%-- Plot figures
for ee = 1:length(inx)
    
    [trials_out{ee}] = ecog_plotSpectra(spectra, gridList(inx{ee}), whichTrials, whichTasks, specs);
    
end

%-- Set axis scale
figlist = get(groot,'Children');
figlist = flipud(figlist(1:length(inx)))';
for ifig=figlist
    axlist = findobj(ifig,'Type','Axes')';
    for iax = axlist
        if ~isempty(iax.Title.String)
            iax.XScale = specs.plot.XScale;
            iax.YScale = specs.plot.YScale;
            if strcmp(iax.XScale,'log')
                if numel(iax.XTick) < 2
                    curtick = [log10(iax.XTick) floor(log10(max(trials.time)))];
                    if diff(curtick) == 0
                        iax.XTick = 10.^([-2:0] + curtick(1));
                    else
                        iax.XTick = 10.^([curtick(1)-diff(curtick) curtick]);
                    end
                elseif numel(iax.XTick) < 3
                    curtick = log10(iax.XTick);
                    iax.XTick = 10.^([curtick(1), mean(curtick), curtick(end)]);
                end
                iax.XTickLabel = cellstr(num2str(iax.XTick',['%.' int2str(log10(diff(iax.XLim))-2) 'f']));
            end
        end
    end  
end

%-- Hide or reorder legend
if isfield(specs,'plot')&&isfield(specs.plot,'showlegend')&&strcmp(specs.plot.showlegend,'no')
    set(findobj(figlist,'Type','Legend'),'Visible','off');  % Hide
else
    for ifig=figlist
        %%% get axis size
        frstaxopos = ifig.Children(end).OuterPosition;
        frstaxpos  = ifig.Children(end).Position;
        
        %%% get index of legend
        legidx = find(strcmp(get(ifig.Children,'Type'),'legend'));
        
        %%% move legend at outer position
        legpos = ifig.Children(legidx).Position;
        ifig.Children(legidx).Position = ...
            [frstaxpos(1)+(frstaxpos(3)-legpos(3))/2, frstaxopos(2)+frstaxopos(4), legpos(3:4)];
        
        %%% move legend and axes next to legend at the top of axes list
        ifig.Children = ifig.Children([legidx:legidx+1,1:legidx-1,legidx+2:end]);
    end
end


