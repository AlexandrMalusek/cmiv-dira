%% Plot results of three-material decomposition
% (Three mass fraction maps)
%
function [] = plotWei3(Wei3, names)

  if nargin < 2
    names{1} = 'Material 1';
    names{2} = 'Material 2';
    names{3} = 'Material 3';
  end
  for i = 1:length(Wei3)
    figure()
    subplot(1,3,1), imagesc(100*Wei3{i}(:, :, 1), [0 100]);
    axis image; axis off; colorbar('SouthOutside'); title(names{1});
    subplot(1,3,2), imagesc(100*Wei3{i}(:, :, 2), [0 100]);
    axis image; axis off; colorbar('SouthOutside'); title(names{2});
    subplot(1,3,3), imagesc(100*Wei3{i}(:, :, 3), [0 100]);
    axis image; axis off; colorbar('SouthOutside'); title(names{3});
  end
end