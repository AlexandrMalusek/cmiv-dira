%% Write the Classified Reconstructed Volume (CRV) to a file
%
% Performs iterative image reconstruction of the Drasim simulated
% anthropomorphic phantom.

clear all
close all
clc

% Options
p2MD = 1;       % using 2MD.
p3MD = 1;       % using 3MS.
numbiter = 4;   % Number of iterations. 

% Function for tissue classification.
tissueClassification = @myTissueClassification;

disp('Setting variables...')
tic

setDiraVariables;

% Run the reconstruction
DIRA
