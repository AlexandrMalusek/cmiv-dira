% Example program that shows how to run the segmentation by the JJ2016 
% algorithm. 
%
% 

%% clear

clearvars; clc;
beep off;  

addpath(genpath('./'))

%% Parameters
dataSet = 2;

atlasFirstSlice = 165;
atlasLastSlice = 195;

%% Select data set
[volume, ImageInfo,sliceThickness] = loadDataSet(dataSet);
scaleInfo = [ImageInfo.PixelSpacing;sliceThickness];
%% Load reference volume
load('/Data/ReferenceImage/ReferenceVolume.mat');
%%
% [ref, ~ ] = loadDataSet(1);
% [ref, refMask] = getReferenceVolume(ref);

%% Load Atlas
atlas = readVisualHuman;
atlas = atlas(:,:,atlasFirstSlice:atlasLastSlice);

atlasVoxelSize = [0.91,0.94,5.0];

%%  %%%%% ------ Segmentation ------ %%%%% %%

[filledCompactBone, ...
 filledMarrow,      ...
 adipose, rectum,   ...
 prostate, air]      = segmentation(volume,scaleInfo,ref,refMask,...
                                    atlas, atlasVoxelSize);
%% Plot function

segmentViewer(volume, filledMarrow, adipose, rectum, prostate, air);

