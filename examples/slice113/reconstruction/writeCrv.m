%% Write the Classified Reconstructed Volume (CRV) to a file
%
% Performs iterative image reconstruction of the Drasim simulated
% anthropomorphic phantom.

clear all
close all
clc

% Function for tissue classification.
tissueClassification = @myTissueClassification;

% Set variables
setDiraVariables;

% Do the reconstruction
DIRA
