function deformedContour = deformableModelProstate(image, atlas)
% deformableModelProstate uses the atlas image to get an initial guess 
% where the prostate is in the CT image. Thereafter the matlab function 
% activecontour is used with the Chan-Vese setting to adapt the initial 
% guess to the image content.
%
% Inputs: CT (gray scale) image of the patient. 
%
%         The atlas: The function segment the prostate from the atlas
%         followed by using the deforable model to find the final position
%         of the prostate. Needs to be either double or single and the
%         values need to be in the range 0 to 255.
% 
% Output: Binary image showing the prostate. 
%

atlas(isnan(atlas))=0;
[~, prostateAtlas]=histc(atlas,[190 210]);
prostateAtlas=bwareaopen(prostateAtlas, 50, 4);

deformedContour = activecontour(image,prostateAtlas,200,'Chan-Vese',6);


end
