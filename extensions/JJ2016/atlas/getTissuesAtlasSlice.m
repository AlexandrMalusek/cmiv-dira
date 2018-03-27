function [ boneTissue, atlasMask, ...
           prostate, rectum]       =   getTissuesAtlasSlice( atlas )
% GETTISSUESATLASSLICE returns binary volumes containing the bones, body,
% prostate and rectum of the atlas "atlasImage". 
%
% Input:    - atlas:        Atlas volume where every slice in the volume 
%                           (z-diection) contains the atlasImage image.
%
% Output:   - bone:         Binary volume containing bone in the atlas 
%                           volume.
%
%           - atlasMask:    Binary volume where all voxels that correspond 
%                           to the body in the atlas are set to one.  
%
%           - prostate:     Binary volume containing the prostate.
%
%           - rectum:       Binary volume containing the rectum. 





% Structuring element
se = ones(3,3,3);

% Create mask containing the body
atlasMask = atlas > 25;

% Get bone tissue
boneTissue = atlas >= 200;
boneRemoveOutlines = imopen(boneTissue,ones(3,3,3)) ;
boneRemoveOutlines = imdilate(boneRemoveOutlines,se,'same');

boneTissue = boneRemoveOutlines & boneTissue;

% Get prostate and rectum
prostate = 125 < atlas & atlas < 175;
prostate = imopen(prostate,se);

rectum   = 75 < atlas & atlas < 125;
rectum   = imopen(rectum,se);



end

