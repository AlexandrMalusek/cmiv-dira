function [deformed_mask]=deformableModelMuscle(image, mask)
% Deformable model used for the segmention of muscles. 
%
% Inputs: Image:    Gray scale image of the patient.
% 
%         mask:  The binary image of the initial position of the muscles.
%
% Output: Binary image containing the segmented muscles. 
%
% EXAMPLE:  deformed_guess = deformableModelMuscle(image, initial_mask);
% 

deformed_mask = activecontour(image,mask,300,'Chan-Vese',7);

end
