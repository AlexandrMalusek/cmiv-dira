function [compactBone, ...
          compactBoneAndBoneMarrow,...
          adipose, ...
          rectum,prostate, ...
          airInsideTheBody, ...
          remainingTissues]      =  segmentation(volume,scaleInfo,...
                                                 ref, refMask,    ...
                                                 atlas, atlasVoxelSize,...
                                                 performRegistration)
                                            
% SEGMENTATION uses the JJ2016 algorithm to segment the bones and the 
% adipose tissue of a CT volume. It also estimates the position of the 
% rectum and the prostate using atlas segmentation.
%
% Inputs:   - volume:              CT volume (double or single).
%
%           - scaleInfo:           Vector containing the voxel size 
%                                  in double.
%
%           - ref:                 Reference volume for the histogram 
%                                  matching. (double or single)
%
%           - refMask:             A binary mask where the voxels 
%                                  corresponding to the body in the 
%                                  reference volume have the value one.
% 
%           - atlas:               The atlas used for the atlas 
%                                  segmentation.
%
%           - atlasVoxelSize:      The voxel size of the atlas (double).
%
%
%           - performRegistration: Variable to turn on and off the 
%                                  registration. Can be set to 'true' or 
%                                  'false'.
%
% 
% The functions outputs are the following binary volumes:
%
%           - compactBone:              The compact bone.
%
%           - compactBoneAndBoneMarrow: the compact bone and the bone 
%                                       marrow.
%
%           - adipose:                  the adipose tissue.
%
%           - rectum:                   the rectum.
%
%           - prostate:                 the prostate.
%
%           - airInsideTheBody:         the air inside the body.

disp('Segmentation ...');
                                                  
[volume, airInsideTheBody, air,bodyMask] = preprocessing(volume);

% Extend the range with values from 0 to 4096 before performing the 
% histogram matching. This is mainly needed when using the JJ2016 algorithm
% together with DIRA.
volumeRemapValues = volume - 10;
volumeRemapValues(volumeRemapValues<0) = 0;
volumeRemapValues = 4096/max(volumeRemapValues(:))*volumeRemapValues;

% Histogram matching
volHistMatch = histogramMatchingWithMask(uint16(volumeRemapValues), ...
                                         uint16(ref), bodyMask, refMask);

                                     
                                     
%------------------------------------%
%    Adipose tissue segmentation     %
%------------------------------------%

[~, adipose] = adiposeSegmentation(volume,air);


%------------------------------------%
%         Bone segmentation          %
%------------------------------------%

% The bones are segmented using threshold segmentation and region growing.
% Not segmented bone tissue is connected by the bone filling algorithm.

% The intervalls in the variable 'thresh' define air and soft tissues; 
% compact bone; compact bone with high image intensity.
thresh = [0 850 1100 1251];   

[bonesThresh,boneSeeds] = thresholdSegmentationBone(volHistMatch,thresh);

boneRegion = regionGrowing(volHistMatch, boneSeeds,  265); % 250

bonesAll = boneRegion | bonesThresh;

[compactBone, compactBoneAndBoneMarrow] = fillingOfBones(bonesAll, volume);



%------------------------------------%
%         Atlas segmentation         %
%------------------------------------%

if(performRegistration)
    disp('Atlas segmentation ...');
    
    bodyMask = bodyMask|airInsideTheBody;

    [resampledVolume,...
     resampledAtlas, ~] = volumeResampling(volHistMatch,atlas, bodyMask,...
                                           scaleInfo,atlasVoxelSize);
                                       
    [ctBoneVolume, ~] = createBoneVolume(resampledVolume);
    
    % When using the VisHum atlas
%     [resampledAtlas, ...
%      boneAtlas, ~] = getTissuesAtlas(resampledAtlas);
    
    % When using atlasImage
    [boneAtlas, ~] = getTissuesAtlasSlice(resampledAtlas);
    

    % Registration algorithm
    registAtlas = registrationMethods(resampledVolume,resampledAtlas,  ...
                                     ctBoneVolume,boneAtlas);
    
    % Get tissues from the VisHum atlas
%     [boneRegion, rectum, prostate] = getTissuesVisHum(registAtlas);

    % Get tissues from the atlasImage
    [boneRegion,~ ,prostate, rectum] = getTissuesAtlasSlice(registAtlas);
    
    % Resample registration results to original scale
    [~, rectum, prostate] = backSamlingResults(boneRegion,rectum, ...
                                              prostate,size(volHistMatch));
    
    rectum   = rectum > 0.7;
    prostate = prostate  > 0.7;
    
else
    rectum   = false(size(volume));
    prostate = false(size(volume));
end




%------------------------------------%
%          Post processing           %
%------------------------------------%

% Resolve tissue conflict if a voxel has been segmented as more than one 
% tissue 
order = {'prostate', 'rectum', 'bone', 'air', 'adipose'};

[compactBoneAndBoneMarrow, ...
 adipose,...
 rectum,prostate, ...
 airInsideTheBody]  =  resolveTissueConflicts(compactBoneAndBoneMarrow, ...
                                              adipose, rectum, prostate,...
                                              airInsideTheBody, order);   

% Find the tissues not segmented by the segmentation algorithms.
segmentedTissues =  compactBoneAndBoneMarrow | adipose | ...
                    rectum | prostate | airInsideTheBody;   
                
remainingTissues =  bodyMask - segmentedTissues;                                               
                                          

end