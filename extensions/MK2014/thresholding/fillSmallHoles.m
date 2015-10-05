function image = fillSmallHoles(original,size)
% Function for filling small holes in a binary image. 
%
% Input: - binary image
%
%        - Integer value corresponding to the max size of the holes. Holes 
%          smaller and equal to the max size are filled by the function.
%       



filled = imfill(original, 'holes');
holes = filled & ~original;
bigholes = bwareaopen(holes, size);
smallholes = holes & ~bigholes;
new = original | smallholes;
image=new;

end