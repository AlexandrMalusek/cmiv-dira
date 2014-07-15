%% Write Drasim generated projections to sinograms

clc;
clear all;
close all;

% Set Matlab path
p = path();
path(p, '../../../functions;../../../data')

%% Output file names
sinogramsFileName = 'sinograms.mat';
sinogramsBhFileName = 'sinogramsBH.mat';
spectraFileName = 'spectra.mat';

% Initialization of variables
M0 = 1152;                      % number of projection angles
N0 = 736;                       % number of detector elements
M1 = 720;	  	        % new nr of projections
N1 = 511;		        % new nr detector elements
dfi1 = 180/M1;		        % new angle increment in degrees
dt0 = 0.067864004196156*pi/180; % detector element arc length [rad]
dt1 = 2*11.525*pi/180/512;      % new detector element size (chosen by Robin)
L  = 0.595;                     % X-ray source to origin distance [m]
dfi0 = 360/M0;                  % angle increment in degrees
gamma = atan(N0/2*dt0/L);       % first angle after rebinning [rad]

writeSinogramsQO;
