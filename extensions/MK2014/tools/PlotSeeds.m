function PlotSeeds(im,seeds,num)
% Function for plotting seed points on an image. 
% 
% Input:    - Gray scale image
%
%           - Seed points is a 2 x N list of the y and the x coordinate of  
%             the seed points position, where N is the number of seeds. 
%
%           - The numbering of the figure
%
% EXAMPLE:
%             seeds = [y_pos1, x_pos1; y_pos2, x_pos2];
%             figure_number = 1;
%
%             PlotSeeds(image,seeds,figure_number)           
%
% Written by Julius Jeuthe, 2015 
figure(num)

% Plot image
imagesc(im); 
colormap gray;
axis image;

% Add seeds
hold on;
for k = 1:size(seeds,1)
    plot(seeds(k,2),seeds(k,1),'b*')
end
hold off;

end