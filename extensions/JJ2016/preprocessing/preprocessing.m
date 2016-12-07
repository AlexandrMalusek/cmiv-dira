function [outputVolume, airInBody, air,bodyMask] = preprocessing(volume)
% Preprocessing of the CT data. All voxels except those of the patient are 
% given the intensity value zero. The function also segments all gray
% values equal or below a certain threshold as air.
%
% Input     - volume:       The CT volume.
%
% Output    - The CT volume. All voxels that are not corresponding to the  
%             patients body are set to zero.
%
%           - Binary volume showing the air inside the body.
%
%           - Binary volume showing the air in the entire CT volume.
%
%           - Binary volume where all voxels corresponding to the body have
%             the value one.

% For DICOM images it is set to 624 (-400 HU)
airThresh = 10.0;

% Threshold segmentation
[xSize,ySize,zSize] = size(volume);
[~, labelVolume]    = histc(volume,[0 airThresh max(volume(:))]); 
labelVolume         = single(labelVolume);

bodyMask = labelVolume == 2;

%% Find biggest object
labeled = bwlabeln(bodyMask,6);

labeledBodyMask   = labeled(bodyMask > 0);
histOfVolSizes    = hist(labeledBodyMask(:),max(labeledBodyMask(:)));
bodyMask          = labeled == (find(max(histOfVolSizes)==histOfVolSizes));


% clean upp
volumeOpened = imopen(bodyMask,ones(5,5,5)); 

%% Fill interier of body
tmpVolume                   = ones(xSize,ySize,zSize+2,'single');
tmpVolume(:,:,2:zSize + 1)  = volumeOpened;
filledBinaryVolume          = imfill(tmpVolume,'holes');

% Removal of the incerted two slices 
filledBodyMask = filledBinaryVolume(:,:,2:zSize+1);

%% Create the output arguments

volume = filledBodyMask.*volume;

air       = volume < airThresh;
airInBody = air & filledBodyMask;
bodyMask  = logical(filledBodyMask - airInBody);

% Convert data to Single
outputVolume = single(volume);

end