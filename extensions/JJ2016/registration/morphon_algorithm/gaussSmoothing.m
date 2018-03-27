function filteredVolume = gaussSmoothing(volume, sigma)
% Function that smooths a volume by a Gaussian filter. 
%
% Input:   - volume:         The 3D volume that is to be smoothed.
%
%          - sigma:          The standard deviation of the Gaussian filter.
%
%
% Output:  - filteredVolume: The smoothed 3D volume.


% Gives an odd filter size
filterSize = round(5*sigma) + mod(ceil(5*sigma),2) - 1;

% Create a Gaussian filter
gaussFilter = fspecial('Gaussian',[filterSize 1],sigma);


% Filter the data
volume = imfilter(volume, gaussFilter,  'replicate');
volume = imfilter(volume, gaussFilter', 'replicate');

filteredVolume = imfilter(volume, permute(gaussFilter, [3 2 1]), ...
                          'replicate');


end