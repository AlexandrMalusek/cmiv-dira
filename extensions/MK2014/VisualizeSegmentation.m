function VisualizeSegmentation(E,I,color)
% Plot the segmentation results on the active image. 
%
% Input: - Image to plot the results on
%        - Segmented binary image
%        - Color of the segmented region. Possible colors are
%          'red', 'green', 'blue', 'yellow', 'margenta' and 'cyan'


 I(find(I~=0))=1;

% Make a truecolor all-green image.
switch color
    case 'red'
        green = cat(3, 2*ones(size(E)), zeros(size(E)), zeros(size(E)));
    case 'green'
        green = cat(3,  zeros(size(E)),0.5*ones(size(E)), zeros(size(E)));
    case 'blue'
        green = cat(3, zeros(size(E)), zeros(size(E)), 2*ones(size(E))); 
    case 'yellow'
        green = cat(3, 0.8*ones(size(E)), 0.8*ones(size(E)), zeros(size(E)));   
    case 'magenta'
        green = cat(3, 1*ones(size(E)), zeros(size(E)), 1*ones(size(E)));
    case 'cyan'
        green = cat(3,  zeros(size(E)),1*ones(size(E)), 1*ones(size(E)));
    otherwise
        green = cat(3, zeros(size(E)), zeros(size(E)), 1*ones(size(E)));
end
hold on
h = imshow(green);
hold off

% Use our influence map image as the AlphaData for the solid
% green image.
set(h, 'AlphaData', I)

end