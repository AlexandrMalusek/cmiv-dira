function maFr = elMaFrImage(pmd, Z, iter)
  % Return elemental mass fractions for element Z directly
  % calculated from mass fractions of doublets and triplets.
  %
  % Input:
  % pmd:  PhantomModelData object
  % Z:    atomic number
  % iter: iteration number (0, ...)
  %
  % Output:
  % maFr: matrix of mass fractions

  nTissueDoublets = length(pmd.tissue2Set{1});
  nTissueTriplets = length(pmd.tissue3Set{1});
  maFr = zeros(size(pmd.recLowSet{1}));

  for it = 1:nTissueDoublets
    for ic = 1:2
      maFr = maFr + pmd.matDoublet{it, ic}.W(Z) * pmd.Wei2Set{iter+1}{it}(:,:,ic);
    end
  end

  for it = 1:nTissueTriplets
    for ic = 1:3
      maFr = maFr + pmd.matTriplet{it, ic}.W(Z) * pmd.Wei3Set{iter+1}{it}(:,:,ic);
    end
  end
end
