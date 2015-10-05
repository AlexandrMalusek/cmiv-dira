function out = findMinMax(binaryImage)
% The function is create a struct that stores the coordinates of the mass 
% point of the binary image as well as the coordinates of the highest and 
% lowest pixel with value 1. 
%
% Input: Binary image
%
% The stored variables in the struct are:
%
%   center_x, center_y:    x and y coordinate of the center
%   max_x, min_x:          highest and lowest x coordinate that's 1
%   max_y, min_y:          highest and lowest y coordinate that's 1
%   
% 
% EXAMPLE:
% struct = findMinMax(binary_image)



[body_y, body_x] = find(binaryImage == 1);  

out.center_x = round(mean(body_x)); 
out.center_y = round(mean(body_y)); 

%% Border

erode_binary = imerode(binaryImage,ones(3,3));
border = binaryImage - erode_binary;

% dist_border = logical(dist_border);
[~,cols_x] = find(border == 1);
[rows_y,~] = find(border == 1);

out.max_x = max(cols_x);
out.min_x = min(cols_x);
out.max_y = max(rows_y);
out.min_y = min(rows_y);

