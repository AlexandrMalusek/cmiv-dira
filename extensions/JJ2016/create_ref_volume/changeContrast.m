function changedVol = changeContrast(vol)
% Function that changes contrast using a linear transformation. The
% parameters were tested on 120 kV images only.  

upperThresh = 1250; % Default: 1250
lowerThresh = 850;  % Default: 850

vol(vol<0) = 0;

changedVol = vol;
mask = changedVol < lowerThresh;
changedVol(mask) = lowerThresh;

mask = changedVol > upperThresh;
changedVol(mask) = upperThresh;

threshDiff = upperThresh-lowerThresh;

changedVol = (upperThresh/threshDiff)*(changedVol - lowerThresh);

end