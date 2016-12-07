function [xPosition,yPosition] = getMassCenter(binaryVolume)
% Function that calculates the mass center of a binary image.
% 
% Example: [xPosition,yPosition] = getMassCenter(binaryVolume);

[y,x,~] =  ind2sub(size(binaryVolume),find(binaryVolume == 1));

xPosition = sum(x(:))/length(x);
yPosition = sum(y(:))/length(y);

end