function bidsEcogPlotElectrodesOnMesh(projectDir, subject, atlas, specs)

% should: 
% - read the electrode positions from the electrodes.tsv. Or should it take
%   an electrode table and surface file as input (read in with
%   bidsEcogReadElectrodeFile and bidsEcogReadSurfFile functions? in which
%   case it should probably be called ecog_plotElectrodesOnMesh (because
%   first input is not projectDir)
%
% - plot an empty mesh with just the electrodes if no atlas specified

% - plot electrodes on top of specified atlas(es) (one figure per atlas)
%   atlas colormaps can be obtained from getAtlasLabels.m

% - plot HDgrid electrodes smaller than regular electrodes

% - have an option to plot the matched node?

% - have an option to colorcode the electrodes according to some value
%   specified as an array of length electrodes?
%
% - have an option to colorcode mesh nodes underlying the electrodes 
%   according to some value and some spread function?

end