%% Write the Classified Reconstructed Volume (CRV) to a file
%
% Performs iterative image reconstruction of the Drasim simulated
% anthropomorphic phantom.

clear all
close all
clc

% Function for tissue classification.
tissueClassification = @myTissueClassification;

rampWindowForMeasuredProjections = @rampWindowForMeasuredProjectionsDefault;
% rampWindowForMeasuredProjections = @myRampWindowForMeasuredProjections;

reconstructIteratedProjections = @reconstructIteratedProjectionsDefault;
% reconstructIteratedProjections = @myReconstructIteratedProjections;

reconstructMeasuredProjections = @reconstructMeasuredProjectionsDefault;
% reconstructMeasuredProjections = @myReconstructMeasuredProjections;

% Set variables
pSetMaterialData = 1;
setDiraVariables;

% Do the reconstruction
DIRA

