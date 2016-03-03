% Run the DIRA post-processing code

setMatlabPath;
load('pmd.mat');
load('smd.mat');
global useCode;
useCode = 0;

% Prostate
prostStr = 'H0.629495C0.117335N0.012196O0.239088Na0.000265P0.000394S0.000571Cl0.000344K0.000312';
prostDens = 1.030;

% Calcium
caStr = 'Ca1.000000';
caDens = 1.550;

% Water
waterStr = 'H0.666667O0.333333';
waterDens = 1.000;

% post-processing triplet
name3SA{1} = 'prostate';
Dens3SA(1) = prostDens;
Att3SA(:, 1) = [Dens3SA(1) * CalculateMAC(prostStr, pmd.eEL),...
  Dens3SA(1) * CalculateMAC(prostStr, pmd.eEH)];

name3SA{2} = 'water';
Dens3SA(2) = waterDens;
Att3SA(:, 2) = [Dens3SA(2) * CalculateMAC(waterStr, pmd.eEL),...
  Dens3SA(2) * CalculateMAC(waterStr, pmd.eEH)];

name3SA{3} = 'calcium';
Dens3SA(3) = caDens;
Att3SA(:, 3) = [Dens3SA(3) * CalculateMAC(caStr, pmd.eEL),...
  Dens3SA(3) * CalculateMAC(caStr, pmd.eEH)];

% Define the prostate mask
load('prostateMask.mat');
maskSA = maskProst;

Wei3SA{1} = MD3(0.01*pmd.recLowSet{5}, 0.01*pmd.recHighSet{5}, Att3SA, Dens3SA, maskSA, 0);
% Plot computed mass fractions from MD3
plotWei3(Wei3SA, name3SA);
WeiAv = MD3SP(mean(0.01*pmd.recLowSet{5}(maskSA)), mean(0.01*pmd.recHighSet{5}(maskSA)), Att3SA, Dens3SA);
fprintf('Average mass fraction m1 = %f, m2 = %f and m3 = %f\n', WeiAv);
