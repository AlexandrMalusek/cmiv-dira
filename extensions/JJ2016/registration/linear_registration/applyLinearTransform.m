function deformedAtlas = applyLinearTransform(atlas, deformation, ...
                                              method, newSize)

% Function that applies a displacement field to the atlas. If the 
% displacemet field was calculated for an atlas that was larger then the 
% volume the atlas was registered to, the size of the volume needs to be 
% given as inargument in the function (see variable 'newSize').
%
% Input:    - atlas:        The atlas volume.
%
%           - deformation:  The displacement field. The 4th dimension  
%                           is used for y, x and z. 
%
%           - method:       Interpolation method. Can be set to 'nearest', 
%                           'linear' or 'cubic'.
%
%           - newSize:      The size of the volume the atlas was registered
%                           to. 
%
%
% Output:   - deformedVol   The deformed atlas.

if nargin == 3
    volSize = size(atlas);
    lengthDiffInZ = 0;
elseif nargin == 4
    volSize = newSize;
    atlasSize = size(atlas);
    lengthDiffInZ = max(atlasSize(3) - volSize(3), 0);
end

[y,x,z] = ndgrid(-floor((volSize(1)-1)/2) : ceil((volSize(1)-1)/2),...       
                 -floor((volSize(2)-1)/2) : ceil((volSize(2)-1)/2),...
                 -floor((volSize(3)-1)/2) : ceil((volSize(3)-1)/2 ...
                 + lengthDiffInZ));



deformedAtlas = ba_interp3(x,y,z,double(atlas), ...
                           x + deformation(:,:,:,2),   ...
                           y + deformation(:,:,:,1),   ...
                           z + deformation(:,:,:,3), method);  
                            
deformedAtlas(isnan(deformedAtlas)) = 0;                            

end