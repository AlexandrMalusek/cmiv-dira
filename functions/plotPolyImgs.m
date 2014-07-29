function [] = plotPolyImgs(recLow, recHigh, iter)
  % plotPolyImgs Plot reconstructed images
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
  axis image;colorbar('horiz');axis off;
  title(sprintf('80kV, recon. No: %d', iter));
  colormap('bone');
  
  subplot(1,2,2);
  imagesc(recHigh, [minHigh maxHigh]);
  axis image;colorbar('horiz');axis off;
  title(sprintf('Sn140kV, recon. No: %d', iter));
  colormap('bone');
end
