function [imageBinaryFilled]= createBinary(image)
% This function creates a binary image (mask) defining the patient's body. 
% The input image is segmented by thresholding followed by hole filling. 
%
% The MK2014 uses this function to get a binary image showing the body of
% the patient in the CT image and the atlas. Thereafter the atlas is
% matched to the patients body by translating and scaling it. 
%
% Input: A grayscale image with format double or single. The function 
%        expect the values to be in the range 0 and 255. 
%
% Output: Binary image showing the patients body.

[~, imageBinary]=histc(image,[7 256]);
imageBinaryFilled=imfill(imageBinary);
end