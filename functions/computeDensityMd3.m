function [densityMd3] = computeDensityMd3(Wei3it, Dens3it)
  % computeMd3Density Compute mass density from the MD3 mass fractions
  %
  % Input:
  % Wei3it:  [Nr x Nr x 3 double] array of mass fractions for a triplet
  % Dens3it: [1 x 3 double] tabulated mass densities of a triplet

  sz = size(Wei3it);
  recRho = zeros(sz(1), sz(2));
  for i = 1:3
    recRho = recRho + Wei3it(:,:,i) / Dens3it(i);
  end
  densityMd3 = 1.0 ./ recRho;
  densityMd3(densityMd3 == Inf) = 0.0;
end
