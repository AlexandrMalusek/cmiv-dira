function [] = plotRecLacImages(recLow, eEL, recHigh, eEH, iter)
  % plotRecLacImages Plot reconstructed images
  %
  % Plot maps of linear attenuation coefficients

  % Scaling for plottning
  minLow = 17;
  maxLow = 33;
  minHigh = 15;
  maxHigh = 23;
  
  figure();
  subplot(1,2,1);
  imagesc(recLow, [minLow maxLow]);
  axis image;colorbar('horiz'); axis off;
  title(sprintf('LAC (1/m), E=%.1f keV, Ni=%d', eEL, iter), 'fontsize', 11);
  colormap('bone');
  
  subplot(1,2,2);
  imagesc(recHigh, [minHigh maxHigh]);
  axis image;colorbar('horiz'); axis off;
  title(sprintf('LAC (1/m), E=%.1f keV, Ni=%d', eEH, iter), 'fontsize', 11);
  colormap('bone');
end
