function newVolume = removeZeroPadding(volume,sizeDiff,volNr)
% Function for removing the zeropadding done by zeroPadVolume.
% 
% Inputs: - volume:     The zero padded volume by zeroPadVolume.
%
%         - sizeDiff:   The difference in size between the input volumes.
%                       This variable is given by the function
%                       zeroPadVolume.
%
%         - volNr:      Value that specifies if the volume was the first or  
%                       the second output from ZeroPadVolume.
%
%   See also zeroPadVolume.


if volNr == 1 
    sizeDiff = -sizeDiff;
end

[sy,sx,sz] = size(volume);

if volNr == 2
    halfDiff = -sizeDiff/2;
    
    startY = 1 +  max(floor(halfDiff(1)),0);
    endY   = sy - max(ceil(halfDiff(1)),0);
    
    startX = 1 +  max(floor(halfDiff(2)),0);
    endX   = sx - max(ceil(halfDiff(2)),0);
    
    startZ = 1;
    endZ   = sz;
      
    newVolume = volume(startY:endY,startX:endX,startZ:endZ);
elseif volNr == 1
    halfDiff = sizeDiff/2;

    startY = 1 +  max(floor(halfDiff(1)),0);
    endY   = sy - max(ceil(halfDiff(1)),0);
    
    startX = 1 +  max(floor(halfDiff(2)),0);
    endX   = sx - max(ceil(halfDiff(2)),0);

    newVolume = volume(startY:endY,startX:endX,:);
end


end