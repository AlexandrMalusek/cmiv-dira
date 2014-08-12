% Print mean values and standard deviations in circular ROIs

setMatlabPath;

% Analyze reconstructed phantom model data
PhantomModelData;
load('pmd.mat');

% Define three circular ROIs
centerX = [148, 272, 335];
centerY = [334, 334, 248];
radius = [20, 20, 20];
color = {'r', 'g', 'b'};

% Print means and standard deviations in ROIs for 4th iteration.
[meanReg, stdReg] = pmd.meanAndStdInRegions(centerX, centerY, radius, color, 4);
disp(meanReg);
disp(stdReg);
