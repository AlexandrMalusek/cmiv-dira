function [isotropVolume, msg] = getIsotropicVoxels(volume, xSize, ...
                                                   ySize,  zSize, method)
% Function that resamples a volume in order for the volume to get isotropic
% voxels. The voxels width in x direction is used as the new voxel size in    
% x, y and z direction. 
% 
% Input:    - volume:    The 3D volume that we wish to resample. 
%
%           - xWidth:    The size of the voxels in x direction.
%
%           - yWidth:    The size of the voxels in y direction.
%
%           - zWidth:    The size of the voxels in z direction.
%
%           - method:    The method used for the resampling. It can be set 
%                        to: 'nearest', 'linear' and 'cubic'.
%
%
% Output:   - isotropVolume:    Volume where there the voxels are resampled
%                               in order to be isotropic. 
%
%           - msg:              Message telling if resampling was needed or
%                               not. If message isn't used as output, the 
%                               function displace a message if no 
%                               resampling was needed. 

[sy, sx, sz] = size(volume);


if((xSize == ySize) && (ySize == zSize))
    % Message
    msg = 'No resampling needed!';
    if nargout < 2
        disp(msg)
    end
    
    isotropVolume = volume;
    return;
    
elseif(xSize == ySize) 
    % Use xWidth as reference size
    zSizeFactor = xSize/zSize;
    
    % Create a grid for the resampling. 
    [y, x, z] = ndgrid(linspace(1,sy,sy),...                    
                       linspace(1,sx,sx),...                    
                       linspace(1,sz,ceil(sz/zSizeFactor)));    

else
    zSizeFactor = xSize/zSize;
    ySizeFactor = xSize/ySize;
    
    [y, x, z] = ndgrid(linspace(1,sy,ceil(sy/ySizeFactor)),... 
                       linspace(1,sx,sx),...                   
                       linspace(1,sz,ceil(sz/zSizeFactor)));       
end

volume = double(volume);

isotropVolume = ba_interp3(volume, x, y, z,method);

isotropVolume = single(isotropVolume);

% Message
msg = 'Volume successful resampled';

end