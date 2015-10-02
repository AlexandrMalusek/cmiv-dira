function [deformed_contour]=deformableModelMuscle(image, contour)
% Deformable model used in order to segment the muscles. For this the
% matlab function activecontour where used with the Chan-Vese setting.
%
% Inputs: CT image   
%         The binary image showing the initial guess of the muscles.

deformed_contour = activecontour(image,contour,300,'Chan-Vese',7);

end
