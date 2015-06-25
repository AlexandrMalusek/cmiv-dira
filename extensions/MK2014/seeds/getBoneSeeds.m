function [seeds] = getBoneSeeds(bonesLabeled)
  % Return seed pixels in the bone

  %create empy seed array
  seeds = zeros(2*max(max(bonesLabeled)), 2);
  
  for i = 1:max(max(bonesLabeled))
    %create binary image from separate parts, one per iteration
    binary = bonesLabeled;
    binary(binary ~= i) = 0;
    binary(binary == i) = 1;
    
    potentialSeeds = bwulterode(binary, 4);  %erode image
    noSeeds = bwmorph(potentialSeeds, 'clean'); %remove single pixels
    
    %check if there are single pixels
    if (max(max(potentialSeeds ~= noSeeds)))
      potentialSeeds = potentialSeeds - noSeeds; %save single pixels
    else
      [row, col] = find(potentialSeeds); %save all remaining pixels
      potentialSeeds = zeros(size(bonesLabeled));
      potentialSeeds(row(1), col(1)) = 1;
    end
    
    %save one seed (the first in list) per separate bone part
    [row,col] = find(potentialSeeds);
    
    if length(row) >= 2
      seeds(i, :) = [row(round(end/3)) col(round(end/3))];
      seeds(max(max(bonesLabeled)) + i, :)= [row(end) col(end)];
    else
      seeds(i, :)=[row(1) col(1)];
    end
    
  %     figure, imagesc(binary), colormap jet;
  %     figure, imagesc(potential_seeds), colormap jet;
  %plot(seeds(:,2),seeds(:,1),'r*');
  end
  
  removeRows = find(seeds == 0);
  if (length(removeRows) > 0)
    for i = 1:(length(removeRows)/2)
      seeds(removeRows(i)-(i-1), :) = [];
    end
  end
end
