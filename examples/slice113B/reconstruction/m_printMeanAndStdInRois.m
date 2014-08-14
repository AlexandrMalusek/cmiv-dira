% Print mean values and standard deviations in circular ROIs

setMatlabPath;
load('pmd.mat');

% Define three circular ROIs
centerX = [266, 113, 358];
centerY = [354, 314, 176];
radius = [20, 20, 20];
color = {'r', 'g', 'b'};

% Print means and standard deviations in ROIs for 4th iteration.
[meanReg, stdReg] = pmd.meanAndStdInRegions(centerX, centerY, radius, color, 4);
disp(meanReg);
disp(stdReg);
