%deformable model
function [deformedContour] = deformableModelProstate(image, atlas)
  % Return a contour of the prostate obtained using the deformable model
	 
  atlas(isnan(atlas)) = 0;
  [~, prostateAtlas] = histc(atlas, [190 210]);
  prostateAtlas = bwareaopen(prostateAtlas, 50, 4);
  
  deformedContour = activecontour(image, prostateAtlas, 200, 'Chan-Vese', 6);
end
