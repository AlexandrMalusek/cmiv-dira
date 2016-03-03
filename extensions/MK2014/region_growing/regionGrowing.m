% The original code by Stéphane was edited by M.K. and A.M.
% Copyright (c) 2012, Stéphane
% All rights reserved.

function Phi = regionGrowing(tolerance, Igray, x, y)
  
  Phi = false(size(Igray, 1), size(Igray, 2));
  PhiOld = Phi;
  Phi(uint16(x), uint16(y)) = 1;
  Phi = imdilate(Phi, strel('disk', 1, 0)); %multipixel seeds
  while(sum(Phi(:)) ~= sum(PhiOld(:)))
    PhiOld = Phi;
    segmVal = Igray(Phi);
    meanSeg = mean(segmVal);
    posVoisinsPhi = imdilate(Phi, strel('disk', 1, 0)) - Phi;
    voisins = find(posVoisinsPhi);
    valeursVoisins = Igray(voisins);
    Phi(voisins(valeursVoisins > meanSeg - tolerance & valeursVoisins < meanSeg + tolerance)) = 1;
  end
end
