%% Plot results of two-material decomposition
% (two mass fraction maps and one mass density map)
%
function [] = plotWei2Dens(Wei2, dens, names)
  
  if nargin < 3
    names{1} = 'Material 1';
    names{2} = 'Material 2';
    names{3} = 'Density';
  end
  for i = 1:length(Wei2)
    figure()
    subplot(1,3,1), imagesc(100*Wei2{i}(:, :, 1), [0 100]);
    axis image; axis off; colorbar('SouthOutside'); title(names{1});
    subplot(1,3,2), imagesc(100*Wei2{i}(:, :, 2), [0 100]);
    axis image; axis off; colorbar('SouthOutside'); title(names{2});
    subplot(1,3,3), imagesc(dens{i});
    axis image; axis off; colorbar('SouthOutside'); title(names{3});
  end
end