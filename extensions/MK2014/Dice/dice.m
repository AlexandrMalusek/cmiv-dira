function [diceValue] = dice(image, groundTruth)
  % Return the Dice similarity coefficient between the input image and the
  % ground truth.

  diceValue = 2*nnz(image & groundTruth) / (nnz(image) + nnz(groundTruth));
end
