%% Calculate polyenergetic projection
%
function [Ap] = computePolyProj(E, uE, N, p, mu)
  
  sizeE = size(E);
  sizeP = size(p);
  
  sl = zeros(sizeP(1), sizeP(2), sizeE(1)-1);
  
  for k = 2:sizeE-1;
    tmpSum = zeros(size(p(:, :, 1)));
    for i = 1:sizeP(3)
      tmpSum = tmpSum+(-mu(E(k), i)*100.*p(:, :, i));
    end
    sl(:, :, k) = (E(k)*N(k))*(E(k+1)-E(k-1)).*exp(tmpSum);
  end
  up = sum(sl, 3)/2;
  
  Ap = -log(up/uE);
end
