function [imageBinaryFilled]= createBinary(image)
%create binary image
[~, imageBinary]=histc(image,[7 256]);
imageBinaryFilled=imfill(imageBinary);
end