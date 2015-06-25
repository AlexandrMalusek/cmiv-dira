function [bones, adipose, prostate, muscles] = segmentation(image, atlas, histogramReference)
  % This is an automated segmentation algorithm.
  % Inputs:
  %     image: the image that will be segmented
  %     atlas: the atlas used for the segmentation
  %     histogramReference: the image used for the histogram matching
  %
  % Outputs:
  %     bones: the bones of the CT-slice represented by a binary image
  %     adipose: the adipose tissue of the CT-slice represented by a binary image
  %     prostate: the prostate of the CT-slice represented by a binary image
  %     muslces: the muscles of the CT-slice represented by a binary image

  %--------------------------------------------------------------------------
  %HISTOGRAM-MATCHING--------------------------------------------------------
  
  %match image to reference image
  [imageHistMatch] = histogramMatching(image, histogramReference);

  %--------------------------------------------------------------------------
  %REMOVE-CT-TABLE-----------------------------------------------------------
  
  %remove CT-table
  imageHistMatch = removeTable(imageHistMatch);
  
  %--------------------------------------------------------------------------
  %THRESHOLDING--------------------------------------------------------------
  
  %thresholding & labeling
  [bonesThreshold, adiposeThreshold, bonesThresholdFilled]...
  = automatedThresholding(imageHistMatch);
  
  %--------------------------------------------------------------------------
  %SEED-GENERATION-----------------------------------------------------------
  
  %get bone seeds
  [boneSeeds] = getBoneSeeds(bonesThreshold);
  
  %get adipose tissue seed
  [adiposeSeed] = getAdiposeSeed(adiposeThreshold);
  
  %--------------------------------------------------------------------------
  %REGION-GROWING------------------------------------------------------------
  
  %region growing on adipose tissue
  adiposeRegionGrowing = regionGrowing(20, imageHistMatch,...
				       adiposeSeed(1), adiposeSeed(2));
  adiposeRegionGrowingFilled = fillSmallHoles(adiposeRegionGrowing,180);
  se = strel('disk',7);
  adiposeRegionGrowingFilled = imclose(adiposeRegionGrowingFilled, se);
  adiposeRegionGrowingFilled = fillSmallHoles(adiposeRegionGrowingFilled, 1000);
  adipose = adiposeRegionGrowingFilled;
  
  %region growing on bones
  bonesRegionGrowing = zeros(size(imageHistMatch));
  for i = 1:length(boneSeeds)
    bonesTemp{i} = regionGrowing(35, imageHistMatch,...
      boneSeeds(i,1), boneSeeds(i,2)); %Bone right hip joint
    bonesRegionGrowing = bonesRegionGrowing | bonesTemp{i};
  end
  bonesRegionGrowingFilled = imfill(bonesRegionGrowing, 'holes');
  
  %--------------------------------------------------------------------------
  %COMBINE TH & RG-----------------------------------------------------------
  
  %bones
  bonesAll = bonesRegionGrowingFilled | bonesThresholdFilled;
  
  %--------------------------------------------------------------------------
  %CLOSE HOLES---------------------------------------------------------------
  
  %closing region growing || thresholding segmentation for bones
  bonesSkeleton = bwmorph(bonesAll, 'skel', Inf);
  se = strel('disk', 9);
  bonesSkeletonClosed = imclose(bonesSkeleton, se);
  bones = bonesSkeletonClosed | bonesAll;
  bones = imfill(bones, 'holes');
  
  %--------------------------------------------------------------------------
  %ATLAS-BASED-IMAGE-REGISTRATION--------------------------------------------
  atlas = rgb2gray(atlas);
  
  warning('off', 'MATLAB:unknownElementsNowStruc');
  [atlasRegistered, atlasRegisteredAffine] =...
  registration(imageHistMatch, atlas, bones);
  warning('on', 'MATLAB:unknownElementsNowStruc');
  
  %--------------------------------------------------------------------------
  %DEFORMABLE MODEL PROSTATE-------------------------------------------------
  
  prostate = deformableModelProstate(imageHistMatch, atlasRegistered);
  
  %--------------------------------------------------------------------------
  %SEGMENTATION OF MUSCLES---------------------------------------------------
  
  [~, atlasMuscles] = histc(atlasRegisteredAffine, [150 180]);
  atlasMuscles = bwareaopen(atlasMuscles, 50, 4);
  atlasMusclesLabeled = bwlabel(atlasMuscles, 4);
  
  muscleLeft = atlasMusclesLabeled;
  muscleLeft(muscleLeft ~= 1) = 0;
  
  muscleRight = atlasMusclesLabeled;
  muscleRight(muscleRight ~= 2) = 0;
  muscleRight(muscleRight == 2) = 1;
  
  muscleLeftFinal = deformableModelMuscle(imageHistMatch.*((uint8(bones)*(-1))+1), muscleLeft);
  muscleRightFinal = deformableModelMuscle(imageHistMatch.*uint8((int8(bones)*(-1))+1), muscleRight);
  
  muscles = muscleLeftFinal + muscleRightFinal;
  
end
