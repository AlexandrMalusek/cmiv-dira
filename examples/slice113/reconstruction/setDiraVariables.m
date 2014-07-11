%% Set DIRA's variables
%
% Set Matlab path
p = path();
path(p, '../../../functions;../../../data')

% Names of files
resultsFileName = 'results.mat';
sinogramsFileName = 'sinograms.mat';
sinogramsBhFileName = 'sinogramsBH.mat';
spectraFileName = 'spectra.mat';
prostateMaskFileName = 'prostateMask.mat';

%% Scanner parameters
%
% X-ray tube voltages in kV
eL = 80;       % low
eH = 140;      % high

% Photon efective energies in keV
eEL = 50.0;    % low
eEH = 88.5;    % high

% CT scanner geometry
L = 0.595;                       % distance source - rot. center in m
alpha = 38.4;                    % fanbeam angle in deg
N1 = 511;                        % number of detector elements after rebinning
dt1 = 0.402298392584693/(N1+1);  % detector element size = pixel distance

% Joseph projection generation
interpolation = 2;

%% Elemental material composition (number of atoms per molecule) and
% mass density (in g/cm^3).

% Compact bone
boneStr = 'H0.403076C0.149395N0.033840O0.316004Na0.001473Mg0.000929P0.034249S0.001056Ca0.059978';
boneDens = 1.920;

% Bone marrow mixture
marrowMixStr = 'H0.620994C0.250615N0.008329O0.119144Na0.000124P0.000092S0.000267Cl0.000240K0.000145Fe0.000051';
marrowMixDens = 1.005;

% Lipid
lipidStr = 'H0.621918C0.341890O0.036192';
lipidDens = 0.920;

% Proteine
proteineStr = 'H0.480981C0.326573N0.089152O0.101004S0.002291';
proteineDens = 1.350;

% Prostate
prostStr = 'H0.629495C0.117335N0.012196O0.239088Na0.000265P0.000394S0.000571Cl0.000344K0.000312';
prostDens = 1.030;

% Calcium
caStr = 'Ca1.000000';
caDens = 1.550;

% Water
waterStr = 'H0.666667O0.333333';
waterDens = 1.000;

% Tissue 1
Dens2{1}(1) = boneDens;
Cross2{1}(:, 1) = [CalculateMAC(boneStr, eEL),...
  CalculateMAC(boneStr, eEH)];
Att2{1}(:, 1) = Dens2{1}(1)*Cross2{1}(:, 1);
mu2Low{1}(:, 1) = CalculateMACs(boneStr, 1:eL);
mu2High{1}(:, 1) = CalculateMACs(boneStr, 1:eH);

Dens2{1}(2) = marrowMixDens;
Cross2{1}(:, 2) = [CalculateMAC(marrowMixStr, eEL),...
  CalculateMAC(marrowMixStr, eEH)];
Att2{1}(:, 2) = Dens2{1}(2)*Cross2{1}(:, 2);
mu2Low{1}(:, 2) = CalculateMACs(marrowMixStr, 1:eL);
mu2High{1}(:, 2) = CalculateMACs(marrowMixStr, 1:eH);

% Tissue 2
Dens3{1}(1) = lipidDens;
Att3{1}(:, 1) = [Dens3{1}(1)*CalculateMAC(lipidStr, eEL),...
  Dens3{1}(1)*CalculateMAC(lipidStr, eEH)];
mu3Low{1}(:, 1) = Dens3{1}(1)*CalculateMACs(lipidStr, 1:eL);
mu3High{1}(:, 1) = Dens3{1}(1)*CalculateMACs(lipidStr, 1:eH);

Dens3{1}(2) = proteineDens;
Att3{1}(:, 2) = [Dens3{1}(2)*CalculateMAC(proteineStr, eEL),...
  Dens3{1}(2)*CalculateMAC(proteineStr, eEH)];
mu3Low{1}(:, 2) = Dens3{1}(2)*CalculateMACs(proteineStr, 1:eL);
mu3High{1}(:, 2) = Dens3{1}(2)*CalculateMACs(proteineStr, 1:eH);

Dens3{1}(3) = waterDens;
Att3{1}(:, 3) = [Dens3{1}(3)*CalculateMAC(waterStr, eEL),...
  Dens3{1}(3)*CalculateMAC(waterStr, eEH)];
mu3Low{1}(:, 3) = Dens3{1}(3)*CalculateMACs(waterStr, 1:eL);
mu3High{1}(:, 3) = Dens3{1}(3)*CalculateMACs(waterStr, 1:eH);

% Tissue 3
Dens3SA(1) = prostDens;
Att3SA(:, 1) = [Dens3SA(1)*CalculateMAC(prostStr, eEL),...
  Dens3SA(1)*CalculateMAC(prostStr, eEH)];
mu3LowSA(:, 1) = Dens3SA(1)*CalculateMACs(prostStr, 1:eL);
mu3HighSA(:, 1) = Dens3SA(1)*CalculateMACs(prostStr, 1:eH);

Dens3SA(2) = waterDens;
Att3SA(:, 2) = [Dens3SA(2)*CalculateMAC(waterStr, eEL),...
  Dens3SA(2)*CalculateMAC(waterStr, eEH)];
mu3LowSA(:, 2) = Dens3SA(2)*CalculateMACs(waterStr, 1:eL);
mu3HighSA(:, 2) = Dens3SA(2)*CalculateMACs(waterStr, 1:eH);

Dens3SA(3) = caDens;
Att3SA(:, 3) = [Dens3SA(3)*CalculateMAC(caStr, eEL),...
  Dens3SA(3)*CalculateMAC(caStr, eEH)];
mu3LowSA(:, 3) = Dens3SA(3)*CalculateMACs(caStr, 1:eL);
mu3HighSA(:, 3) = Dens3SA(3)*CalculateMACs(caStr, 1:eH);

% Ordering of coefficients are important first tissues for 2MD in the
% output order of tissue classification and then tissues for 3MD in output
% order of tissue classification.
tissueOrder2 = [1];
tissueOrder3 = [1];

muLow = cat(2, mu2Low{1}, mu3Low{1});
muHigh = cat(2, mu2High{1}, mu3High{1});
