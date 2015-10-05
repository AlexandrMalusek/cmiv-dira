function diceValue = dice(image, groundTruth)
% The function dice calculates the dice cooficient. A dice cooficient of 1
% corresponds to a perfect match and a dice cooficient of 0 cooresponds to
% that no pixels in the binary images matches the ground truth.
%
% Input: Segmentation result. For instance the binary image of the bones
%
%        The ground thruth needs is a binary image showing the true 
%        position of the tissue.  
%
% Output: Dice cooficient. 

diceValue = 2*nnz(image & groundTruth)/(nnz(image) + nnz(groundTruth));

end