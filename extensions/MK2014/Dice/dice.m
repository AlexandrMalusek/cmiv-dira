function diceValue = dice(image, groundTruth)
% The function dice calculates the dice coefficient. A dice coefficient of 1
% corresponds to a perfect match and a dice coefficient of 0 cooresponds to
% no pixels in the binary image matching the ground truth.
%
% Input: Segmentation result. For instance the binary image of the bones
%
%        The ground thruth needs is a binary image showing the true 
%        position of the tissue.  
%
% Output: Dice coefficient. 

diceValue = 2*nnz(image & groundTruth)/(nnz(image) + nnz(groundTruth));

end