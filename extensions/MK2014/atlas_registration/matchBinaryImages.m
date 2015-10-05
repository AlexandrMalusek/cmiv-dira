function atlasMoved = matchBinaryImages(BinaryImage, BinaryAtlas,atlas)
% Match two binary images by moving the binary region of the atlas to the 
% position of the patient body. The binary regions are matched by first 
% calculating the lowest and highest value of the binary region (body of 
% the patient) in x and y direction followed by translating and scaling the
% atlas to match those.
% 
%
% Inputs:   - Binary image: The binary representation of the patients body. 
%
%           - Binary atlas: The binary representation of the body in the 
%                    atlas.
%
%           - Atlas: The atlas that should be moved by the transformation.
%                    It is expected that the atlas is of type double or 
%                    single. It is also expected that all three inputs
%                    have the same image size.
%
% Output:   - Moved Atlas: Atlas matched according to the algorithm.


% Create a structure element with the min and max position of the images 
im_str = findMinMax(BinaryImage);
at_str = findMinMax(BinaryAtlas);

% Get the scaling factor in x and y direction
ratio_x = (at_str.max_x - at_str.min_x)/(im_str.max_x - im_str.min_x); 
ratio_y = (at_str.max_y - at_str.min_y)/(im_str.max_y - im_str.min_y); 

% Translate and rescale the image. In order to make the scaling easy the
% left corner is used as the reference point for the scaling. The binary
% image is moved there, rescaled and thereafter placed on the correct
% position
new_atlas = atlas(at_str.min_y:end,at_str.min_x:end);

new_atlas = padarray(new_atlas,[at_str.min_y,at_str.min_x],'post');

[nsy,nsx] = size(new_atlas);
[x,y]   = meshgrid(1:nsx,1:nsy);

[x2,y2]   = meshgrid((1 : ratio_x :ratio_x*(nsx+1)), ...
                     (1 : ratio_y : ratio_y*(nsy+1)));  

% interpolate to the new size
atlasMoved = interp2(x,y,double(new_atlas),x2,y2);
atlasMoved(isnan(atlasMoved))=0;

% Give the binart image the right position
[sy,sx] = size(BinaryAtlas);
atlasMoved = padarray(atlasMoved,[im_str.min_y,im_str.min_x],'pre');
atlasMoved = atlasMoved(1:sy,1:sx);


end