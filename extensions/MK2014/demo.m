% Plot a CT image segmented using the MK2014 algorithm and the corresponding
% ground truth. Print the Dice similarity coefficient. The ground truth was
% created by en experienced radiologist using manual segmentation. The MK2014
% algorithm needs an atlas image (atlasImage_CT01.png) for the atlas-based
% segmentation and a reference image (referenceImage_CT01.jpg) for the
% histogram matching. They have to be provided by the user.

close all;
clear all;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
tic
%recursively add sub-folders to path
addpath(genpath('./'));

%testimages
image = double(dicomread('image_1'));

%atlas
atlas = imread('atlasImage_CT01.png');

%histogram reference image
ref = imread('referenceImage_CT01.jpg');

[bones, adipose, prostate, muscles] = segmentation(image, atlas, ref);

load('groundTruthSegmentation_CT01.mat')
diceBones = dice(bones, bonesGT);
diceProstate = dice(prostate, prostateGT);
diceAdipose = dice(adipose, adiposeGT);
diceMuscles = dice(muscles, musclesGT);

X = ['Dice Values:',...
    ' Bones: ', num2str(diceBones),...
    ' Prostate: ', num2str(diceProstate),...
    ' Adipose: ', num2str(diceAdipose),...
    ' Muscles: ', num2str(diceMuscles)];
disp(X)
%--------------------------------------------------------------------------
%PLOTS---------------------------------------------------------------------

%segmented image
figure, imagesc(image), colormap gray;
VisualizeSegmentation(image, muscles, 'blue');
VisualizeSegmentation(image, prostate, 'red');
VisualizeSegmentation(image, bones, 'green');
VisualizeSegmentation(image, adipose, 'yellow');

%ground truth
figure, imagesc(image), colormap gray;
VisualizeSegmentation(image, musclesGT, 'blue');
VisualizeSegmentation(image, prostateGT, 'red');
VisualizeSegmentation(image, bonesGT, 'green');
VisualizeSegmentation(image, adiposeGT, 'yellow');
toc
