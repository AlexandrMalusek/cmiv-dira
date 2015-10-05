function [deformed_contour]=deformableModelMuscle(image, contour)
% Deformable model used in order to segment the muscles. For this the
% matlab function activecontour where used with the Chan-Vese setting.
%
% Inputs: Gray scale image showing the patient.
% 
%         The binary image showing the initial guess of the muscles.
%
% Output: Binary image showing the deformed initial guess. 
%
% EXAMPLE:  deformed_guess = deformableModelMuscle(image, initial_guess);
% 

deformed_contour = activecontour(image,contour,300,'Chan-Vese',7);

end
