%% Set DIRA's variables
%

setMatlabPath;

% Define what code to use
% 0 = Matlab, 1 = C, 2 = OpenMP
global useCode
useCode = 0;

%% Load data and initialize variables
% ------------------------------------
disp('Loading data and initializing variables...')
tic

% Names of files
sinogramsFileName = 'sinograms.mat';
sinogramsBhFileName = 'sinogramsBH.mat';
spectraFileName = 'spectra.mat';

%% Scanner model data
smd = ScannerModelData;
smd.eL = 80;              % low x-ray tube voltage in kV
smd.eH = 140;             % high x-ray tube voltage in kV
smd.L = 0.595;            % distance source - rot. center in m
smd.N1 = 511;             % number of detector elements after rebinning
smd.dt1 = 0.402298392584693/(smd.N1+1); % detector element size
smd.interpolation = 2;    % Joseph projection generation
smd.N0 = 512;             % nr of detector elements
smd.M0 = 560;             % nr of projections
smd.fact = 4;             % factor for no rebinned projections
smd.dfi0 = 228 / smd.M0;  % angle increment in degrees
smd.dfi1 = 1 / smd.fact;  % new angle increment in degrees
smd.M1 = 180 * smd.fact;  % new nr of projections
smd.dt0 = (38.4 * pi / 180) / smd.N0;  % Detector element arc length [rad]
smd.gamma = atan(smd.N0 / 2 * smd.dt0 / smd.L);  % First angle after rebinning [rad]

% Phantom model data
pmd = PhantomModelData;
pmd.savedIter = [0, 1, 2, 4]; % save data for these iterations
pmd.eEL = 50.0;       % low effective energy in keV
pmd.eEH = 88.5;       % high effective energy in keV
pmd.p2MD = 1;         % using 2MD.
pmd.p3MD = 1;         % using 3MS.
pmd.recAlg = 0;       % reconstruction algorithm (0 = old DIRA, 1 = iterative DEFBP)

% The setting of material data takes long time. Skip it if not needed.
if ~pSetMaterialData
  return
end

spectra =  load(spectraFileName);
smd.ELow = spectra.currSpectLow(1:75, 1);
smd.NLow = spectra.currSpectLow(1:75, 2);
smd.EHigh = spectra.currSpectHigh(1:135, 1);
smd.NHigh = spectra.currSpectHigh(1:135, 2);

sinograms = load(sinogramsFileName);
pmd.projLow = sinograms.projLow;
pmd.projHigh = sinograms.projHigh;

sinogramsBH = load(sinogramsBhFileName);
pmd.projLowBH = sinogramsBH.projLowBH;
pmd.projHighBH = sinogramsBH.projHighBH;

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

% Protein
proteinStr = 'H0.480981C0.326573N0.089152O0.101004S0.002291';
proteinDens = 1.350;

% Water
waterStr = 'H0.666667O0.333333';
waterDens = 1.000;

% Material doublets
pmd.matDoublet{1,1} = Material('compact bone', boneDens, boneStr);
pmd.matDoublet{1,2} = Material('bone marrow', marrowMixDens, marrowMixStr);

% Material triplets
pmd.matTriplet{1,1} = Material('lipid', lipidDens, lipidStr);
pmd.matTriplet{1,2} = Material('protein', proteinDens, proteinStr);
pmd.matTriplet{1,3} = Material('water', waterDens, waterStr);

