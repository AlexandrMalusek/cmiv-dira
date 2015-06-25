function [imageBinaryFilled] = createBinary(image)
  % Create binary image

  [~, imageBinary] = histc(image, [7 256]);
  imageBinaryFilled = imfill(imageBinary);
end
