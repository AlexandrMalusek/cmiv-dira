%% Two-material decomposition in DECT.
%
function [Wei2, dens] = MD2(AttE1mat, AttE2mat, Att2, Dens2, mask)
  %
  % Input:
  % AttE1mat: matrix of measured LACs at effective energy E1
  % AttE2mat: matrix of measured LACs at effective energy E2
  % Att2:     tabulated LACs of doublet base materials at E1 and E2
  % Dens2:    mass density of doublet base materials
  % mask:     mask defining the tissue to be decomposed
  %
  % Output:
  % Wei2:     2 matrices of mass fractions
  % dens:     matrix of mass densities
  %

  imgSize = size(mask);
  Wei2 = zeros(imgSize(1), imgSize(2), 3);
  dens = zeros(imgSize(1), imgSize(2));
  
  for k = 1:imgSize(1)
    for l = 1:imgSize(2)
      if mask(k, l) == 1
        M = [(Att2(1, 1)/Dens2(1)-Att2(1, 2)/Dens2(2)) -AttE1mat(k, l);
          (Att2(2, 1)/Dens2(1)-Att2(2, 2)/Dens2(2)) -AttE2mat(k, l)];
        b = -[Att2(1, 2)/Dens2(2); Att2(2, 2)/Dens2(2)];
        w = M\b;
        
        Wei2(k, l, 1) = w(1);
        Wei2(k, l, 2) = 1-w(1);
        
        dens(k, l) = 1/w(2);
      end
    end
  end
end
