%Thresholding
function [bonesLabeled, adiposeLabeled, bonesFilledNoNoise]...
  = automatedThresholding(image)

  %threshold image at set ranges
  [~, imageThreshold] = histc(image, [0 6 75 180 256]);
  
  %--------------------------------------------------------------------------
  %BONES---------------------------------------------------------------------
  %separate bones
  bones = imageThreshold;
  bones(bones ~= 4) = 0;
  %remove small objects from the bone image
  bonesNoNoise = double(bwareaopen(bones, 50, 4));
  %fill bones
  bonesFilledNoNoise = fillSmallHoles(bonesNoNoise, 10000);
  
  %label each different bone part as an individual integer
  bonesLabeled = bwlabel(bonesNoNoise, 4);
  
  %--------------------------------------------------------------------------
  %ADIPOSE-------------------------------------------------------------------
  
  %separate adipose tissue
  adipose = imageThreshold;
  adipose(adipose ~= 2) = 0;
  %remove small objects from adipose image
  adiposeNoNoise = double(bwareaopen(adipose, 50, 4));
  
  %label each different adipose tissue part as an individual integer
  adiposeLabeled = bwlabel(adiposeNoNoise, 4);

  %remove all adipose tissue except the largest area
  adiposeLabeledHistogram = imhist(uint8(adiposeLabeled));
  adiposeLabeledHistogram(1) = 0;
  [~, largestArea] = max(adiposeLabeledHistogram);
  largestArea = largestArea - 1;
  adiposeLabeled(adiposeLabeled ~= largestArea) = 0;
  
end
