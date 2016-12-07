function [vol, atlas,sizeDiff] = zeroPadVolume(vol,atlas,sizeDiff)

% Function for zero padding the volume (vol) and the atlas to the same size
% (voxel array dimensionality) in x and y direction. When the volume has
% more elements in z-direction than the atlas, the volume is also zero 
% padded in z-direction.
%   
% Input:    - vol:      The volume
%
%           - atlas:    The atlas
%
%           - sizeDiff: Vector containing the size difference between the 
%                       volumes (optional). If the size difference is not 
%                       given as input argument, it is calculated by the
%                       function.  
%
% Output:   - vol:      The zero padded volume
%
%           - atlas:    The zero padded atlas
%
%           - sizeDiff: Vector describing the size difference before the 
%                       zero padding. This vector is needed when removing
%                       the added zeros.
%
% Example: 
% % Zero pad the volumes
% [volPad, atlasPad,sizeDiff]  = zeroPadVolume(vol,atlas);
% [volMaskPad, atlasMaskPad,~] = zeroPadVolume(volMask,atlasMask,sizeDiff);
%
% % Remove the added zeros
% volNr = 1;  % vol was the first input of zeroPadVolume
% vol   = removeZeroPadding(vol,sizeDiff,volNr);
%
% 
% see also removeZeroPadding.

if nargin == 2
%if numel(sizeDiff) == 1 && sum(sizeDiff.^2) == 0
    fixedSize  =  size(vol);
    atlasSize  =  size(atlas);
    
    sizeDiff = fixedSize - atlasSize;
end

% zero-padding of x, y and z
if(sizeDiff(1) < 0)
    halfDiffY = -sizeDiff(1)/2;
    vol = padarray( vol ,floor(halfDiffY),'pre');
    vol = padarray( vol ,ceil(halfDiffY),'post');
elseif(sizeDiff(1) > 0)
    halfDiffY = sizeDiff(1)/2;
    atlas = padarray( atlas ,floor(halfDiffY),'pre');
    atlas = padarray( atlas ,ceil(halfDiffY),'post');
end

if(sizeDiff(2) < 0)
    halfDiffX = -sizeDiff(2)/2;
    vol = padarray( vol ,[0 floor(halfDiffX)],'pre');
    vol = padarray( vol ,[0 ceil(halfDiffX)],'post');
elseif(sizeDiff(2) > 0)
    halfDiffX = sizeDiff(2)/2;
    atlas = padarray( atlas ,[0 floor(halfDiffX)],'pre');
    atlas = padarray( atlas ,[0 ceil(halfDiffX)] ,'post');
end

if(sizeDiff(3) > 0)
    halfDiffZ = sizeDiff(3)/2;
    atlas = padarray( atlas ,[0 0 floor(halfDiffZ)],'pre');
    atlas = padarray( atlas ,[0 0 ceil(halfDiffZ)],'post');
end

end