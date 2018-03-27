function [boneVolume, mask] = createBoneVolume(volume)
% Function that creates a volume that only contains the bone tissue in the 
% CT volume.  
%
% Input:    - volume:       The CT volume.
%
%
% Output:   - boneVolume:   Volume containing the bone tissue. All other
%                           tissues are set to zero.
%
%           - mask:         Binary volume showing the bone tissue. 

boneThresh = 700;%820;

%% Create bone volume
mask = volume < boneThresh;
boneVolume         = volume;
boneVolume(mask)   = boneThresh; 
boneVolume         = boneVolume - boneThresh;

mask = ~mask;

end