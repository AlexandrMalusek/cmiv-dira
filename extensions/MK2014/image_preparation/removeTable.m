function [imageNoTable] = removeTable(image)
  % Remove CT table from the image

  [~, imageBinary] = histc(image, [7 256]);
  imageBinaryFilled = imfill(imageBinary);
  imageBinaryNoTable = double(bwareaopen(imageBinaryFilled, 20000));
  imageNoTable = uint8(imageBinaryNoTable) .* image;
end
