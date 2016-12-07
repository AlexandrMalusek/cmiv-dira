function volumeWithFilledHoles = fillSmallHoles(volume,maxSize,conn)
% Function that fills small holes in a binary volume. A hole is defined as
% a binary region with value zeros surrounded by voxels with value one.
%
% Inputs:   - volume    Binary volume that is to be filled.
%
%           - maxSize   The largest size a hole is allowed to have to be
%                       filled by the algorithm.
%
%           - conn      The connectivity used for the closing. 'conn' can
%                       be set to 4 and 8 for closing in 2D and 6, 18 and
%                       26 for closing in 3D.
%
%
%
% Outputs:  - volumeWithFilledHoles     Binary volume where all holes
%                                       smaller than maxSize have been
%                                       closed.
%
%
% See bwareaopen for more information about the closing.


volumeWithFilledHoles = ~double(bwareaopen(~volume, maxSize, conn));

end
