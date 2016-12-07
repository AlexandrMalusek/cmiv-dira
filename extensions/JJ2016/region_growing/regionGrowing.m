function finalRegion = regionGrowing(volume, seedVolume, meanTol)
% Region growing algorithm that adds a voxel to the seed region if the 
% voxel is ajacent to the seed region and has a gray value inside the 
% acceptance range. The acceptance range is defined as the intervall 
% [regionMean - meanTolerance, regionMean + meanTolerance], where
% regionMean is the mean value of all voxels inside the grown region (seed
% region) and meanTol is a user defined tolerance.
%
%
% Inputs:       - volume:     Volume containing the data set. 
%
%               - seedVolume: Binary volume containing the seed points or
%                             seed region.
%
%               - meanTol:    Value on how much the gray value of an voxel 
%                             is allowed to deviate from the mean value of 
%                             the seed region in order to be included into   
%                             the seed region.
%
% Output:       - finalRegion: The region growing result as a binary
%                              volume.


% Create a Structure Element
structureElement = zeros(3,3,3);
structureElement(1:3,2,2) = 1;
structureElement(2,1:3,2) = 1;
structureElement(2,2,1:3) = 1;

newInterestRegion = seedVolume;
interestRegion    = zeros(size(volume));
lastIteration     = 0;

while(sum(sum(sum(newInterestRegion - lastIteration)))>0)
    % Calculate mean of the new voxels
    regionMean = mean(volume(logical(newInterestRegion)));
    lastIteration = newInterestRegion;

    
    while(sum(interestRegion(:)) ~= sum(newInterestRegion(:)) || stopFlag)
        stopFlag = false;
        interestRegion = logical(newInterestRegion);
        
        % Get border
        dilatedRegion = imdilate(interestRegion,structureElement);
        borderRegion  = logical(dilatedRegion - interestRegion);
        
        acceptedMeanVox = (volume(borderRegion) > regionMean - meanTol)&...
                          (volume(borderRegion) < regionMean + meanTol); 
        % Voxels accepted 
        acceptedVox = acceptedMeanVox;
        
        newInterestRegion(borderRegion == 1) = acceptedVox;  
    end
    stopFlag = true; 
    
end

finalRegion = newInterestRegion;

end

