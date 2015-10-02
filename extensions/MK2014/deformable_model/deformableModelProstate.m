function [deformedContour]=deformableModelProstate(image, atlas)
% Deformable model used in order to segment the prostate. For this the
% matlab function activecontour where used with the Chan-Vese setting.
%
% Inputs: CT image   
%
%         The atlas: The function segment the prostate from the atlas
%         followed by using the deforable model to find the final position
%         of the prostate
%         

atlas(isnan(atlas))=0;
[~, prostateAtlas]=histc(atlas,[190 210]);
prostateAtlas=bwareaopen(prostateAtlas, 50, 4);

deformedContour = activecontour(image,prostateAtlas,200,'Chan-Vese',6);


end
