function visualizeSegmentation(E,seg_result,color,alpha)
% Plot the segmentation results on the active image. 
%
% Input: - Image to plot the results on
%        - Segmented binary image
%        - Color of the segmented region. Possible colors are
%          'red', 'green', 'blue', 'yellow', 'margenta' and 'cyan'

influence_map = alpha*double(seg_result);
influence_map(influence_map<0) = 0;
% influence_map(seg_result~=0) = alpha;

% Make a truecolor all-green image.
switch color
    case 'red'
        true_color = cat(3, ones(size(E)), zeros(size(E)), zeros(size(E)));
    case 'green'
        true_color = cat(3,  zeros(size(E)),0.5*ones(size(E)), zeros(size(E)));
    case 'blue'
        true_color = cat(3, zeros(size(E)), zeros(size(E)), 1*ones(size(E))); 
    case 'yellow'
        true_color = cat(3, 0.8*ones(size(E)), 0.8*ones(size(E)), zeros(size(E)));   
    case 'magenta'
        true_color = cat(3, 1*ones(size(E)), zeros(size(E)), 1*ones(size(E)));
    case 'cyan'
        true_color = cat(3,  zeros(size(E)),1*ones(size(E)), 1*ones(size(E)));
    otherwise
        true_color = cat(3, zeros(size(E)), zeros(size(E)), 1*ones(size(E)));
end
hold on
h = imagesc(true_color); %axis image; axis off
hold off

% Use our influence map image as the AlphaData for the solid
% green image.
set(h, 'AlphaData', influence_map)

end