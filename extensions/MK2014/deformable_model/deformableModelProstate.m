function deformedmask = deformableModelProstate(image, atlas)
% deformableModelProstate uses the atlas image to get an initial position 
% of the prostate. Thereafter the matlab function activecontour is used 
% with the Chan-Vese setting to adapt the initial position to the image 
% content.
%
% Inputs: CT (gray scale) image of the patient. 
%
%         The atlas: An atlas defining the position of the prostate. The 
%         image needs to be in either double or single precision and the 
%         values need to be in the range 0 to 255.
% 
% Output: Binary image containing the prostate. 
%

image(image < 0) = 0;
atlas(isnan(atlas))=0;

[~, prostateAtlas]=histc(atlas,[190 210]);
prostateAtlas=bwareaopen(prostateAtlas, 50, 4);

deformedmask = activecontour(image,prostateAtlas,200,'Chan-Vese',6);


end
