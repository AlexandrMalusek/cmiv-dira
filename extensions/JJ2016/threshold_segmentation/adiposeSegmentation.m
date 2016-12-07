function [adiposeThresh, adiposeFilled] = adiposeSegmentation(volume,air)
% Function that segments the adipose tissue using threshold segmentation. 
% The segmentation is done by a lower fixed threshold and an upper 
% threshold calculated using Otsu's method. The segmentation result is 
% given by the first output argument.
%
% The second output contains a post-processed version of the segmentation 
% result, where voxels that are adjacent to voxels containing air are 
% excluded for the segmentation result. Additionally, holes up to 20 
% voxels are closed in the segmentation result (binary volume). 
%
% Input:    - volume    CT volume
%
%           - air       Binary volume containg the air
%
% Output:
%           - adiposeThresh     Binary volume contaning the adipose tissue.
%
%           - adiposeFilled     Binary volume contaning the the post-
%                               processed segmentation result of the 
%                               adipose tissue. 



% Parameters
threshLow = 10.0;                           % DICOM images: 624 (-400 HU) 

% Highest and lowest intensity value for finding the correct 
% threshold for Otsus method.
threshMin = 18;                            % DICOM images: 824  (-200 HU) 
threshMax = 26;                            % DICOM images: 1124 ( 100 HU) 




% Finding the threshold that separate the adipose tissue from the remaining 
% tissue by using Otsu's method
threshOtsu = otsu(volume,threshMin,threshMax);

adiposeThresh = threshLow< volume & volume<threshOtsu;

% Remove partial volume effect
dialateAir = imdilate(air,ones(3,3,3));

adiposeWithoutPartial = (adiposeThresh - (dialateAir & adiposeThresh)) > 0;


% Filling of holes in the segmentation result. 
adiposeFilled = fillSmallHoles(adiposeWithoutPartial,20,6);

end