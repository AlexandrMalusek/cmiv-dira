function [newBoneMask, filledBoneMarrow]= fillingOfBones(boneMask, ctVolume)
% Algorithm for connecting segments of compact bone and fill the interior  
% of the compact bone. Only voxels that have a higher local image intensity
% than the voxels in their neighborhood are used for connecting the 
% segments of compact bone. Such a connection must additional be connected 
% to the compact bone in at least two different places.
%
% The filling of the interior (segmentation of the bone marrow) is 
% performed by morphological filling, filling each slice separately. The 
% interior of non-filled slices are closed by dilation of the binary volume
% in z-direction. 
%
% Input:    - boneMask:  Binary volume containing the compact bone.
%
%           - ctVolume:  Volume containing the CT data.
%
%
%
% Output:   - newBoneMask:      Binary volume containing the filled compact 
%                               bone.
%
%           - filledBoneMarrow: Binary volume containing the compact bone  
%                               and the bone marrow.


% High pass & high intensity filltering of the CT data to get the local
% maximum mask (LMM).
averagedVol      = smooth3(ctVolume,'box',[9 9 9]); %9
localMaximumMask = ctVolume > averagedVol;

%% Finds possible connections by connecting segments of compact bone
se = ones(3,3,3);
% Dilatation where the growing of the bone mask is restricted by the LMM.   
dilatedMask = double(bwareaopen(boneMask, 50, 4));
for k = 1:5
    dilatedMask = localMaximumMask.*imdilate(dilatedMask,se); 
end

% Define the region where the bones are allowed to be closed.
connectBone = boneMask;

for k = 1:7
    connectBone = imdilate(connectBone,se,'same');
end
for k = 1:7
    connectBone = imerode(connectBone,se,'same');
end

filledCompactBone = dilatedMask .* connectBone;

% Remove small objects in 2D. This makes it easier to remove objects that
% only are connected in the z-direction
removeSmallregions = double(bwareaopen(filledCompactBone, 50, 4));
removeSmallregions = localMaximumMask.*imclose(removeSmallregions,se);
compactBoneWithoutSmallRegions = removeSmallregions | boneMask;


%% Find regions where the filled region is connected to the original 
% compact bone

filledRegions        = compactBoneWithoutSmallRegions - boneMask;
labeledFilledRegions = bwlabeln(filledRegions,26);

% structure element
se2 = zeros(3,3,3); 
se2(1:3,2,2) = 1;
se2(2,1:3,2) = 1; 
se2(2,2,1:3) = 1;

dilatedBoneMask = imdilate(boneMask,se2);

%% Remove all connections
% Remove all connections not connected to two or more different regions of
% compact bone

newBoneMask = boneMask;
for k = 1:max(labeledFilledRegions(:))
    holeRegion    = labeledFilledRegions == k;
    connectRegion = dilatedBoneMask & holeRegion;
    labeledConnectRegion = bwlabeln(connectRegion,26);
    % Histogram is used to determine if the there is an overlap in more
    % than one place
    n = hist(labeledConnectRegion(connectRegion(:)));
    n = n > 0;
    if(sum(n)>1) 
        newBoneMask = newBoneMask | holeRegion;
    end
end


%% Filling the bone marrow in 2D
filledBoneMarrow = zeros(size(newBoneMask));

for k = 1:size(boneMask,3)
   filledBoneMarrow(:,:,k)=imfill(compactBoneWithoutSmallRegions(:,:,k),...
                                  4,'holes');
end

%% Closing in z-direction is used to connect the CT slices that couldn't be
% filled by the slice-wise filling of the interior of the compact bone.

zSize = 7;
filledBoneMarrow = logical(closingZDirection(filledBoneMarrow,zSize));
filledBoneMarrow = filledBoneMarrow | localMaximumMask.* ...
                                    imclose(filledBoneMarrow,ones(5,5,5));

removeSmallregions = double(bwareaopen(filledBoneMarrow, 50, 4));
filledBoneMarrow   = boneMask | removeSmallregions;



end