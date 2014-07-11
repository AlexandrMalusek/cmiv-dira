%% Three-material decomposition in DECT
%
function [Wei3] = MD3(AttE1mat, AttE2mat, Att3, Dens3, mask)

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