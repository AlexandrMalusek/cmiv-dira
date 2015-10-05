function ImageMatched = histogramMatching(image,ref)
% Function for matching the histogram of an image to the histogram of the 
% reference image. 
%
% Input:  - Image: The image on which the histogram matching should be
%           performed on. The image needs to be uint8 or double in the
%           range 0 to 255. 
%
%         - Reference image: The reference image should have the format
%           uint8.
%
% Output: - Histogram matched image. The image is in the format uint8
%
% EXAMPLE: 
%    matched_image = histogramMatching(image,reference_image)

if(size(ref,3)==3)
    ref=rgb2gray(ref);
end
image=uint8(255*(image/(max(max(image)))));

refBig=imresize(ref,4,'bicubic');
imageBig=imresize(image,4,'bicubic');
imageMatchedBig=imhistmatch(imageBig,refBig,256);
ImageMatched=imresize(imageMatchedBig,0.25,'bicubic');

end