function [bone, boneSeed] = thresholdSegmentationBone(volume,thresh)
% Threshold segmentation of the bone tissue using two different thresholds. 
% The lower threshold is used to segment the bone tissue. The second higher
% threshold is used to segment regions of compact bone with high image 
% intenstity. 
%
% Inputs:  - volume:    CT volume
%
%          - thresh:    Vector containing two thresholds for segmentation 
%                       of bone tissue. The first threshold needs to be 
%                       lower than the second one. 
%
% Output:  - bone:      Binary volume containing the voxels above the lower
%                       threshold. 
% 
%          - boneSeed:  Binary volume containing the voxels above the upper
%                       threshold. 
%


[~, imageThreshold]=histc(volume,thresh);

bone     = (imageThreshold == 2) | (imageThreshold == 3);
boneSeed =  imageThreshold == 3;


end