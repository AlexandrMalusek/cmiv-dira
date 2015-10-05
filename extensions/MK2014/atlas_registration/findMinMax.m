function out = findMinMax(binaryImage)
% The function creates a structure that stores the coordinates of the 
% body center in the binary image as well as the largest and the lowest 
% coordinates of the body.
%
% Input: Binary image
%
% The stored variables in the struct are:
%
%   center_x, center_y:    x and y coordinate of the center
%   max_x, min_x:          largest and lowest x coordinate of the body
%   max_y, min_y:          largest and lowest y coordinate of the body
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

