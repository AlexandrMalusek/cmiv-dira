function closedVolume = closingZDirection( vol,filterSizeZ )
% Function for performing closing in z direction. 
%
% Input:    - Binary volume that is to be closed.
%           - The filter size in z-direction
%
% Output:   - The closed binary volume


zBorder =  ceil(filterSizeZ/2);

extendedVolume = padarray(vol,[0 0 zBorder],'both');

extendedVolumeClosed = imclose(extendedVolume,ones(1,1,filterSizeZ)); 
closedVolume = extendedVolumeClosed(:,:,zBorder+1:end-zBorder);


end

