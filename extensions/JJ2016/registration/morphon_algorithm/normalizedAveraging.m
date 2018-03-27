function filteredVolume = normalizedAveraging(volume,filterSize,sigma,certainty) 
% Function that performs normalized averaging of a volume. 
%
% Inputs:   - volume:       The 3D volume that is to be smoothed.
%
%           - filterSize:   The filter size of the Gaussian kernel.
%
%           - sigma:        The standard deviation of the Gaussian kernel.
%
%           - certainty:    The certainty.
%
%
% Output  - filteredVolume: The smoothed volume.



% To avoid division by zero
certainty = certainty + eps;

gaussKernel = fspecial('Gaussian',[filterSize 1],sigma);

% Perform smoothing of the volume.
volume = convn(volume.*certainty,gaussKernel, 'same');
volume = convn(volume, gaussKernel', 'same');
volume = convn(volume, permute(gaussKernel, [3 2 1]), 'same');

% Perform smoothing of the certainty
certainty = convn(certainty, gaussKernel, 'same');
certainty = convn(certainty, gaussKernel', 'same');
certainty = convn(certainty, permute(gaussKernel, [3 2 1]), 'same');

% Set voxels that have a certainty of zero to 1 in order to avoid division
% by 0.
zeroCert = find(certainty==0);
certainty(zeroCert) = 1;
volume(zeroCert) = 0;

% Remove the certainty that was used as weight function
filteredVolume = volume./certainty;


end