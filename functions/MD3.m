%% Three-material decomposition in DECT
%
function [Wei3] = MD3(AttE1mat, AttE2mat, Att3, Dens3, mask)
  %
  % Input:
  % AttE1mat: matrix of measured LACs at effective energy E1
  % AttE2mat: matrix of measured LACs at effective energy E2
  % Att3:     tabulated LACs of triplet base materials at E1 and E2
  % Dens3:    mass density of triplet base materials
  % mask:     mask defining the tissue to be decomposed
  %
  % Output:
  % Wei3:     3 matrices of mass fractions
  %
  % Examples:
  %
  % In the function call:
  % Wei3{i} = MD3(AttE1mat, AttE2mat, Att3{tissueOrder3(i)},...
  %    Dens3{tissueOrder3(i)}, tissue3{i})
  % the following variables are used:
  %
  % >> size(AttE1mat)
  % ans =
  %    511   511
  %
  % >> size(AttE2mat)
  % ans =
  %    511   511
  %
  % >> Att3{tissueOrder3(1)} % LACs (1/cm) of 3 triplet tissues at 2 energies
  % ans =
  %     0.1909    0.2814    0.2270
  %     0.1602    0.2271    0.1772
  %
  % >> Dens3{tissueOrder3(1)} % density (g/cm^3) of 3 triplet tissues
  % ans =
  %     0.9200    1.3500    1.0000
  %
  % >> size(Wei3{1}) % mass fractions of 3 triplet tissues
  % ans =
  %    511   511     3

  imgSize = size(mask);
  Wei3 = zeros(imgSize(1), imgSize(2), 3);
  
  for k = 1:imgSize(1)
    for l = 1:imgSize(2)
      if mask(k, l) == 1
        M = [ (AttE1mat(k, l) - Att3(1, 1))/Dens3(1) - ...
          (AttE1mat(k, l) - Att3(1, 3)) / Dens3(3), ...
          (AttE1mat(k, l) - Att3(1, 2)) / Dens3(2) - ...
          (AttE1mat(k, l) - Att3(1, 3)) / Dens3(3); ...
          (AttE2mat(k, l) - Att3(2, 1)) / Dens3(1) - ...
          (AttE2mat(k, l) - Att3(2, 3)) / Dens3(3), ...
          (AttE2mat(k, l) - Att3(2, 2)) / Dens3(2) - ...
          (AttE2mat(k, l) - Att3(2, 3)) / Dens3(3)  ];
        b = (- [(AttE1mat(k, l) - Att3(1, 3)) / Dens3(3); ...
          (AttE2mat(k, l) - Att3(2, 3)) / Dens3(3)]);
        w = M \ b;
        
        Wei3(k, l, 1) = w(1);
        Wei3(k, l, 2) = w(2);
        Wei3(k, l, 3) = 1 - w(1) - w(2);
      end
    end
  end
end
