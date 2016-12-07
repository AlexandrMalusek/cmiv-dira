function newGrayValues = histogramMatchingWithMask(volume,refVolume,...
                                                   bodyMask, refMask)
% Histogram matching function where the histogram matched regions are
% defined by binary volumes. The voxels not transformed by the histogram
% matching are set to zero.
%
% Inputs:   - volume    The volume for which the gray values are 
%                       transformed.
%
%           - refVolume The reference volume
%
%           - bodyMask  Binary volume defining the region in which the gray
%                       values should be transformed by the histogram 
%                       matching.
%
%           - refMask   Binary volume defining the region used for the 
%                       histogram matching of the reference volume.
%
%
%
% Output:   - newGrayValues     The histogram matched volume.
%


% Set the number of possible gray values, n, dependent on the type of the 
% image data
if(isa(volume,'uint8') && isa(refVolume,'uint8'))
    n = 256;
elseif(isa(volume,'uint16') && isa(refVolume,'uint16'))
    n = 4096; % 65535
    % set to 4096 since most medical images only uses 12 bytes
    % Will need to much memory if we use 65535 steps
else
    error('Input images needs to be uint8 or uint16');
end

% Create histograms
numberOfBodyVoxels = numel(volume(bodyMask(:)));

% Get reference Histogram 
refHistogram = histc(double(refVolume(refMask(:))),0:(n-1)); 

% Normalize the reference histogram and calculate the CDF
refHistogram    = (numberOfBodyVoxels/sum(refHistogram)) * refHistogram;       
cumReferenceVol = cumsum(refHistogram);
lenghtRef       = length(refHistogram);

% Create the histogram and the CDF of the CT volume
histOfVolume = histc(double(volume(bodyMask(:))),0:(n-1));
CDFOfVolume  = cumsum(histOfVolume);


%% Find best match

% Test condition
diff = (cumReferenceVol(:)*ones(1,n) - ones(lenghtRef,1)*CDFOfVolume(:)'); 

% Remove all values that don't fullfil the condition from the evaluation
% by setting them to a high number
mask = diff + numberOfBodyVoxels * sqrt(eps) < 0; 
diff(mask) = numberOfBodyVoxels;

[~,T] = min(diff); 
T = T-1;


%% Transform the values of the volume

if(isa(volume,'uint8'))
    newGrayValues = intlut(volume,uint8(T));
elseif(isa(volume,'uint16'))
    lut = repmat(T,[1 16]);                     % 65536/4096 = 16
    newGrayValues = intlut(volume,uint16(lut(1:65536)));
end

% Ensures that the air in b is zero. In all test cases have shown that this
% isn't necessary but it's alway good to be safe 
newGrayValues(bodyMask == 0) = 0; 


end






