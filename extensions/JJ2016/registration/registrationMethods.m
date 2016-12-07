function registAtlas = registrationMethods(ctVol,atlas,bone,atlasBone)
% Function that performs the atlas segmentation in the JJ2016 algorithm. It
% performs first a linear registration using phase-based affine
% registration followed by non-linear registration by the Morphon. 
%
%
% Input:    - ctVol:     The CT volume.
% 
%           - atlas:     The atlas volume, that is to be registered to the  
%                        CT volume. 
%
%           - bone:      3D volume containing the bone of the CT volume.
%
%           - atlasBone: Binary volume containing the bone of the atlas 
%                        volume.
%
%
% Output:   - registAtlas  The deformed atlas volume.  


% Settings

affineIterations   = [10 10 30 30 60 100];
morphonIterations  = [5, 15, 15, 15, 10, 10, 10, 10];
scalesAffine       = [4 3]; 
scalesMorphon      = [5 3];

% Variable that can be one ore zero. If one, the difference between the
% mass centers of the bone in the CT volume and the atlas is calculated and
% the translation needed to set them to the same mass center set as initial
% transformation for the linear registration.
doInitialtrans = 0;

%% Make volumes the same size

%[bone, atlas_bone,sizeDiff] =  zeroPadVolume(bone,atlas_bone,0);
[bone, atlasBone,sizeDiff] =  zeroPadVolume(bone,atlasBone);


[ctVol, atlas] =  zeroPadVolume(ctVol,atlas,sizeDiff);

%% Create weight function

affineWeight = atlasWeightFunction(bone,floor(size(bone,3)/4),'z');
affineWeight = affineWeight .* atlasWeightFunction(bone,20,'x');
affineWeight = affineWeight .* atlasWeightFunction(bone,20,'y');

%% Set initial translation

if(doInitialtrans == 1)
    compactBone = ctVol > 1100;
    [centerVolX, centerVolY]    = getMassCenter(compactBone);
    [centerAtlasX,centerAtlasY] = getMassCenter(atlasBone); 
    initialTranslation = [(centerVolX   - centerAtlasX) ...
                          (centerAtlasY - centerVolY)];
else
    initialTranslation = [0 0];
end
%% Affine registration

[linearRegisteredAtlas, ...
 displaceField, ~] = linearRegistration(bone, atlasBone, scalesAffine, ...
                                        affineIterations,affineWeight, ...
                                        initialTranslation);
                                    
atlasDeformedByLinear = applyLinearTransform(atlas, ...
                                             displaceField,size(bone));

%% Morphon

[sy,sx,sz] =  size(bone);

croppedBoneAtlas = linearRegisteredAtlas(1:sy,1:sx,1:sz);
croppedAtlas     = atlasDeformedByLinear(1:sy,1:sx,1:sz);


% Use Run the Morphon
[~, accumulatedDisplacement] = morphon(bone, croppedBoneAtlas,  ...
                                       scalesMorphon, morphonIterations);      
                                       
nonLinearDeformedAtlas = applyMorphonDisplacement(croppedAtlas, ...
                                         accumulatedDisplacement);


%% Remove zero padding

registAtlas = removeZeroPadding(nonLinearDeformedAtlas,sizeDiff,2);

end