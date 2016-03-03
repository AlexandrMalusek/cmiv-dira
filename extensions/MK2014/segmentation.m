function [bones,adipose,prostate,muscles] = segmentation(image,atlas,histogramReference)
%This is an automated segmentation algorithm.
%Inputs:
%     image:               - The image that are to be segmented
%
%     atlas:               - The atlas used for the segmentation
%
%     histogramReference:  - The image used for the histogram matching
%
% Outputs:
%     bones:    - The bones of the CT-slice represented by a binary image
%
%     adipose:  - The adipose tissue of the CT-slice represented by a 
%                 binary image
%
%     prostate: - The prostate of the CT-slice represented by a binary 
%                 image
%
%     muslces:  - The muscles of the CT-slice represented by a binary image


% Preprocessing of the image
image(image<500) = 24;

%--------------------------------------------------------------------------
%HISTOGRAM-MATCHING--------------------------------------------------------

%match image to reference image
[imageHistMatch]=histogramMatching(image,histogramReference);
% -------------------------------------------------------------------------
%REMOVE-CT-TABLE-----------------------------------------------------------

%remove CT-table
imageHistMatch=removeTable(imageHistMatch); 

%--------------------------------------------------------------------------
%THRESHOLDING--------------------------------------------------------------
%thresholding & labeling
[bonesThreshold, bonesThresholdFilled,adiposeCleaning]...
    =automatedThresholding(imageHistMatch);

%--------------------------------------------------------------------------
% SEED-GENERATION----------------------------------------------------------

%get bone seeds
[boneSeeds]=getBoneSeeds(bonesThreshold);

%get adipose tissue seed
[adiposeSeed]=getAdiposeSeed(adiposeCleaning);

%--------------------------------------------------------------------------
% REGION-GROWING-----------------------------------------------------------

%region growing on adipose tissue
adiposeRegionGrowing = zeros(size(image));
adiposeRegion = cell(size(adiposeSeed,1),1);

for k = 1:size(adiposeSeed,1)
    adiposeRegion{k} = regionGrowing(20,imageHistMatch,... 
                           adiposeSeed(k,1),adiposeSeed(k,2));
    adiposeRegionGrowing = adiposeRegionGrowing | adiposeRegion{k};
end
adiposeRegionGrowingFilled=fillSmallHoles(adiposeRegionGrowing,180);


se = strel('disk',7);
adiposeRegionGrowingFilled=imclose(adiposeRegionGrowingFilled,se);
adiposeRegionGrowingFilled=fillSmallHoles(adiposeRegionGrowingFilled,1000);
adipose=adiposeRegionGrowingFilled;

%% region growing on bones
bonesRegionGrowing=zeros(size(imageHistMatch));
bonesTemp = cell(size(boneSeeds,1),1);

for i=1:size(boneSeeds,1)
    bonesTemp{i}=regionGrowing(32,imageHistMatch,...
        boneSeeds(i,1),boneSeeds(i,2)); %Bone right hip joint
    bonesRegionGrowing=bonesRegionGrowing | bonesTemp{i};
end
bonesRegionGrowingFilled=imfill(bonesRegionGrowing,'holes');


% -------------------------------------------------------------------------
%COMBINE TH & RG-----------------------------------------------------------

%bones
bonesAll=bonesRegionGrowingFilled | bonesThresholdFilled;

%--------------------------------------------------------------------------
%CLOSE HOLES---------------------------------------------------------------

%closing region growing || thresholding segmentation for bones
bonesSkeleton = bwmorph(bonesAll,'skel',Inf);
se = strel('disk',9);
bonesSkeletonClosed = imclose(bonesSkeleton,se);
bones=bonesSkeletonClosed | bonesAll;
bones=imfill(bones,'holes');


% -------------------------------------------------------------------------
%ATLAS-BASED-IMAGE-REGISTRATION--------------------------------------------

% Ensure that the atlas is an gray scale image  
if size(atlas,3) == 3
    atlas=rgb2gray(atlas);
end

warning('off', 'MATLAB:unknownElementsNowStruc');
% Registrate outlines
[sy, sx] = size(image);
[x,y] = meshgrid(-(sx-1)/2:(sx-1)/2,-(sy-1)/2:(sy-1)/2);
%create binary image and atlas
[imageBinary]= createBinary(imageHistMatch);
[atlasBinary]= createBinary(atlas);

% Simple Registration of Outline
atlasMovedOutline = matchBinaryImages(imageBinary, atlasBinary, atlas);

% Bone registration
iter = [200,60,40,30,30,40,40];
[~, bonesAtlas]=histc(atlasMovedOutline,[220 256]);
[vx,vy] = registration(bones, bonesAtlas, 3,iter);

bonesMovedAffine = interp2(x,y,double(atlasMovedOutline),x+vx,y+vy);
warning('on', 'MATLAB:unknownElementsNowStruc');

% -------------------------------------------------------------------------
% DEFORMABLE MODEL PROSTATE------------------------------------------------

prostate = deformableModelProstate(imageHistMatch, bonesMovedAffine);

% -------------------------------------------------------------------------
% SEGMENTATION OF MUSCLES--------------------------------------------------

[~, atlasMuscles]=histc(atlasMovedOutline,[150 180]);
atlasMuscles=bwareaopen(atlasMuscles, 50, 4);
atlasMusclesLabeled=bwlabel(atlasMuscles,4);

muscleLeft=atlasMusclesLabeled;
muscleLeft(muscleLeft~=1)=0;

muscleRight=atlasMusclesLabeled;
muscleRight(muscleRight~=2)=0;
muscleRight(muscleRight==2)=1;

muscleLeftFinal=deformableModelMuscle(imageHistMatch.*((double(bones)*(-1))+1), muscleLeft);
muscleRightFinal=deformableModelMuscle(imageHistMatch.*double((double(bones)*(-1))+1), muscleRight);

muscles=muscleLeftFinal+muscleRightFinal;


end
