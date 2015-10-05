function [imageBinaryFilled]= createBinary(image)
% This function creates a binary image showing the outline of the body. The
% body is segmented by thresholding followed by hole filling. 
%
% The MK2014 uses this function to get a binary image showing the body of
% the patient in the CT image and the atlas. Thereafter the atlas is
% matched to the patients body by translating and scaling it. 
%
% Input: A grayscale image. The thresholds in the functions are set to 
%        match the histogram matched image.
%
% Output: Binary image showing the patients body.

[~, imageBinary]=histc(image,[7 256]);
imageBinaryFilled=imfill(imageBinary);
end