function [backsampledVolume, msg] = backSampling(vol,oldVolumeSize, method)
% Function for resampling the volume to its initial size. It can be used to 
% resample a volume with isotropic voxels to its original scale and should 
% be used together with the function 'getIsotropicVoxels'.
%
% Input:    - vol:           3D volume that needs to be resampled.
%
%           - oldVolumeSize: The old volume size of the 
%
%           - method:        The method used for the interpolation. The
%                            options are 'nearest', 'linear' and 'cubic'.
%
% Output:   - backsampledVolume:   The resampled volume.
%
%           - msg:                 Message containing information about if
%                                  the volume needed to be resampled or 
%                                  not. 
%
%
%   See also getIsotropicVoxels.



[sy, sx, sz] = size(vol);

if(size(vol) == oldVolumeSize)
    msg = 'No resampling needed for this volume!';
    if nargout < 2
        disp(msg)
    end
    
    backsampledVolume = vol;
    return;
    
elseif(sy == sx)
   
    % create a coordinate system
    [y, x, z]=...
                ndgrid(linspace(1,sy,sy),...               
                       linspace(1,sx,sx),...              
                       linspace(1,sz,oldVolumeSize(3)));  
    
else   
    % create a coordinate system
    [y, x, z]=...
                ndgrid(linspace(1,sy,oldVolumeSize(1)),...
                       linspace(1,sx,sx),...              
                       linspace(1,sz,oldVolumeSize(3)));
    
end

backsampledVolume = interp3(vol, x, y, z, method);

msg = 'The volume was successful resampled to its original size!';

end