function [imageNoTable]=removeTable(image)
% Removes all objects but the biggest object (the patient) in the image.
% The area surrounding the patient is considered air and set to -100 to
% increase the contrast between adispose tissue and air.
%
% Input:  CT image: The image must either have the format double or single 
%                   and have values in the range 0 to 255. 
%                   
% Output: CT image without table and air set to -100. Input format is the
%         same as input


%remove CT-table
[~, imageBinary]=histc(image,[7 256]);
imageBinaryFilled=imfill(imageBinary);

% Removal of all obejcts but the biggest one
imageBinaryNoTable=bwlabel(imageBinaryFilled,4);
imageLabeledHistogram=imhist(uint8(imageBinaryNoTable));
imageLabeledHistogram(1)=0;
[~, largestArea]=max(imageLabeledHistogram);
largestArea=largestArea-1;
imageBinaryNoTable(imageBinaryNoTable~=largestArea)=0;

% imageBinaryNoTable = double(bwareaopen(imageBinaryFilled, 20000));
imageNoTable=imageBinaryNoTable.*double(image);

imageNoTable(imageNoTable==0) = -100;

end