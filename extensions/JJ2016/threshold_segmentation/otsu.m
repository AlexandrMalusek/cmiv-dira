function thresh = otsu(volume,lowerThresh,upperThresh)
% Functions that finds a threshold using Otsu's method of a 3D volume. Only
% voxels with a gray value in between the lower and the upper threshold
% are used for the calculations.
%
% Inputs:   - volume:       The volume for which the threshold is 
%                           calculated. 
%
%           - lowerThresh:  Lower threshold value.
%
%           - upperThresh:  Upper threshold value.
%
%
%
% Output:   - thresh:       The threshold value given by Otsu's method.   
%




% Gray values in DIRA are not limited to be integer  
steps = 100;

% Find voxels 
newVolume    = volume >= lowerThresh & volume <= upperThresh;
nrTotalVoxel = sum(newVolume(:));
volume2      = newVolume.*volume;

counts = hist(volume2(:),steps*ceil(max(volume2(:))));

% Rescale the thresholds for DIRA.
upperThresh = steps*upperThresh;
lowerThresh = steps*lowerThresh;

% Initialization for loop
sumWeightBack = 0;
maxBetweenVar = 0;


for iter = lowerThresh : upperThresh-1
    % Calculate the weights
    sumWeightBack = sumWeightBack + counts(iter);
    sumWeightFor  = nrTotalVoxel - sumWeightBack;
    
    weightBack = sumWeightBack/nrTotalVoxel;
    weightFor  = sumWeightFor/nrTotalVoxel;
    
    
    % Mean for background and forground
    meanBack = sum((lowerThresh:iter) .* counts(lowerThresh:iter)) ...
               / sumWeightBack;
           
    meanFor  = sum((min(iter+1,upperThresh):upperThresh).* ...
             counts(min(iter+1,upperThresh):upperThresh))/sumWeightFor;
    
         
    % Between class variance     
    betweenClassVariance = weightBack*weightFor*(meanBack - meanFor).^2;
   
    % Check if a new maxima is found
    if (betweenClassVariance > maxBetweenVar)
        maxBetweenVar = betweenClassVariance;
        thresh = iter;
    end
end

thresh = thresh / steps;

end