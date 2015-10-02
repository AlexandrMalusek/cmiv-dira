function [ImageMatched]=histogramMatching(image,ref)
% Matches the histogram of the image to the histogram of the reference 
% image. 
%
% EXAMPLE: 
%    matched_image = histogramMatching(image,reference_image)

ref=rgb2gray(ref);
image=uint8(255*(image/(max(max(image)))));

refBig=imresize(ref,4,'bicubic');
imageBig=imresize(image,4,'bicubic');
imageMatchedBig=imhistmatch(imageBig,refBig,256);
ImageMatched=imresize(imageMatchedBig,0.25,'bicubic');

end