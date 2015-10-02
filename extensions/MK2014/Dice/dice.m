function [diceValue]=dice(image, groundTruth)
% Calculates the dice cooficient of the image. 
%
% Input: Segmentation result. For instance the binary image of the bones
%
%        The Ground Thruth for the tissue, e.g. bones
%
% Output: Dice value
diceValue = 2*nnz(image & groundTruth)/(nnz(image) + nnz(groundTruth));

end