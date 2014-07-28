%% Write Drasim generated projections to sinograms

clc;
clear all;
close all;

% Set Matlab path
p = path();
path(p, '../reconstruction;../../../functions;../../../data')

%% Output file names
sinogramsFileName = 'sinograms.mat';
sinogramsBhFileName = 'sinogramsBH.mat';
spectraFileName = 'spectra.mat';

pSetMaterialData = 0;  % Material data are not needed.
setDiraVariables;

writeSinograms;
