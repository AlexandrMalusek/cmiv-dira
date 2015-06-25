function [seed] = getAdiposeSeed(adipose)
  % Return a seed pixel in the adipose tissue

  se = strel('disk', 3);
  adipose = imerode(adipose, se);
  potentialSeeds = bwulterode(adipose, 4);  %erode image
  
  noSeeds = bwmorph(potentialSeeds, 'clean'); %remove single pixels
  
  %check if there are single pixels
  if (max(max(potentialSeeds ~= noSeeds)))
    diff = potentialSeeds - noSeeds; %save single pixels
  else
    [row, col] = find(potentialSeeds); %save all remaining pixels
    diff = zeros(size(bones_labeled));
    diff(row(1), col(1)) = 1;
  end

  [row, col] = find(diff);
  if row(1) == 1, row(1) = 3; end
  if col(1) == 1, col(1) = 3; end
  
  seed=[row(1) col(1)];
end
