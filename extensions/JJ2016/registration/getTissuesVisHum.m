function [boneTissue, rectum, prostate] = getTissuesVisHum(deformedAtlas)
% Function that finds the bones, the rectum and the prostate in the 
% deformed VisHum atlas. The function returns the three tissues as binary 
% volumes, where the tissues are set to one. 
% 
% Input:    - deformedAtlas:    The deformed atlas volume.
%
%
% Output:   - boneTissue:       The bone tissue.
%
%           - rectum:           The rectum.
%
%           - prostate:         The prostate.


% Structuring element
se = ones(3,3,3);

% Get bone segmentation from registration
boneTissue = deformedAtlas >57 & deformedAtlas <118;
boneRemoveOutlines = imopen(boneTissue,ones(3,3,3)) ;
boneRemoveOutlines = imdilate(boneRemoveOutlines,se,'same');

boneTissue = boneRemoveOutlines & boneTissue;

% Get prostate and rectum
prostate = round(deformedAtlas) == 49;
prostate = imopen(prostate,se);

rectum   = round(deformedAtlas) == 23;
rectum   = imopen(rectum,se);


end