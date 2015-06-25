function [deformed_contour] = deformableModelMuscle(image, contour)
  % Return a contour of the muscle obtained using the deformable model

  deformed_contour = activecontour(image, contour, 300, 'Chan-Vese', 7);
end
