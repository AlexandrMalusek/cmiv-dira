function [labeledBones, bonesFilled,adiposeRegion]...
    =automatedThresholding(image)
% Function for performing automated thresholding of a CT image. It is
% assumed that the CT image is enhanced by the histogram matching. If
% another reference image is used than the image that is default in the
% MK2014, the threshold values will need to be modified. 
% 
% Outputs: - labeledBones:   Is used in order to get the seeds of the bones 
%
%          - bonesFilled:    Binary image showing the position of the bones
%
%          - adiposeRegion:  Binary adipose image. Is later used to get the
%                            seeds for the region growing

%threshold image at set ranges
[~, imageThreshold]=histc(image,[0 6 75 180 256]);

%--------------------------------------------------------------------------
%BONES---------------------------------------------------------------------
%separate bones
bones=imageThreshold;
bones(bones~=4)=0;
%remove small objects from the bone image
bonesNoNoise = double(bwareaopen(bones, 50, 4));
%fill bones
bonesFilled=fillSmallHoles(bonesNoNoise,10000);

% Using a higher threshold for finding the seeds
[~, compactBoneThreshold]=histc(image,[0 225 256]);

bones=compactBoneThreshold;
bones(bones~=2)=0;

bonesNoNoise = double(bwareaopen(bones, 50, 4));

%label each different bone part as an individual integer
labeledBones = bwlabel(bonesNoNoise,4);


% --------------------------------------------------------------------------
%ADIPOSE-------------------------------------------------------------------

%separate adipose tissue
adipose=imageThreshold;
adipose(adipose~=2)=0;
%remove small objects from adipose image
adiposeNoNoise = double(bwareaopen(adipose,50,4));

% label each different adipose tissue part as an individual integer
labeledAdipose=bwlabel(adiposeNoNoise,4);

%remove all adipose tissue except the largest area
adiposeLabeledHistogram=imhist(uint8(labeledAdipose));
adiposeLabeledHistogram(1)=0;
[~, largestArea]=max(adiposeLabeledHistogram);
largestArea=largestArea-1;
labeledAdipose(labeledAdipose~=largestArea)=0;

% Do some cleaning of the binary adipose image. By eroding the image we
% ensure that the seeds aren't on the edges of the tissue. 
adiposeRegion = double(bwareaopen(adipose,200,4));
adiposeRegion = imerode(adiposeRegion,ones(3,3));

end


