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

%% Sinograms
N0      = 512;               % nr of detector elements
M0      = 560;               % nr of projections

%% Geometry
L = 0.595;                   % distance source - rot. center
fact = 4;                    % factor for no rebinned projections
dfi0 = 228/M0;               % angle increment in degrees
dt1 = 2*11.525*pi/180/N0;    % new detector element size
dfi1 = 1/fact;               % new angle increment in degrees
N1 = 511;                    % new nr detector elements
M1 = 180*fact;               % new nr of projections

% Settings of drasim
%-------------------
alpha = 38.4;                % Fan beam angel [degrees]
gLen = alpha*pi/180;         % Detector arc lenght [rad]
dt0 = gLen/N0;               % Detector element arc length [rad]
gamma = atan(N0/2*dt0/L);    % First angle after rebinning [rad]

writeSinograms;
