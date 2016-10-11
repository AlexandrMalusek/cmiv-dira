function [seed]=getAdiposeSeed(adipose)
% Function for getting a vector with seed points from a binary image. The
% binary regions are first eroded once to remove spurs followed by
% an erosion to point. 
%
% For the picking of the seeds the image is divided to 16 equally sized
% parts. In each part one seed is taken.
%
% Input:  - Binary image showing the adipose tissue
%
% Output: - The output is a 2 x N vector with the y and x coordinate of the 
%           seed points, where N is the number of found seed points.

blockSize = 128;

se = strel('disk',3);
adipose=imerode(adipose,se);
potentialSeeds=bwulterode(adipose, 4);  %erode image

noSeeds = bwmorph(potentialSeeds,'clean'); %remove single pixels

imSize = size(adipose);
numberBlocks = ceil(imSize/blockSize);

seed = [];
for kx = 1:numberBlocks(2)
    for ky = 1:numberBlocks(1)
        intx = (1+blockSize*(kx-1)):min((blockSize*kx),imSize(2));
        inty = (1+blockSize*(ky-1)):min((blockSize*ky),imSize(1));
        
        potSeeds = potentialSeeds(inty,intx);
        noSeed = noSeeds(inty,intx);
        
        %check if there are single pixels
        if(max(max(potSeeds~=noSeed)))
            diff=potSeeds - noSeed; %save single pixels
        else
            diff = potSeeds; 
        end
        
        
        [row,col]=find(diff);
        if(~isempty(col))
            seed = [seed;(128*(ky-1)+row(1)) (128*(kx-1)+col(1))];
        end
    end
end


end