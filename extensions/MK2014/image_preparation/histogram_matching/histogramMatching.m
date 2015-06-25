function [ImageMatched] = histogramMatching(image, ref)
  % Return an image whose intensities are changed so that the histograms of
  % intensities of this and the reference image match.

  ref = rgb2gray(ref);
  image = uint8(255*(image/(max(max(image)))));
  
  refBig = imresize(ref, 4, 'bicubic');
  imageBig = imresize(image, 4, 'bicubic');
  imageMatchedBig = imhistmatch(imageBig, refBig, 256);
  ImageMatched = imresize(imageMatchedBig, 0.25, 'bicubic');
end
