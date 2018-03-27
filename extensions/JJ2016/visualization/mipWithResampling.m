function im = mipWithResampling(vol, dir, resamp, voxelSize, ...
                                plotImage, figNumber)
% Projects the maximum values of a volume on a plane. 
%
% Inputs: 
%   - volume        The volume.
%
%   - direction:    Can be set to 'x', 'y' or 'z'
%
%   - resamp:       Resamples the volume to isotropic voxels if set to 1.
%                   
%   - voxelSize:    Vector containing the size of the voxel in x, y and z
%                   direction.
%
%   - plotImage:    If set to 1, the function plots the calculated MIP 
%                   image.
%
%   - figNumber:    The figure number of the image. This is only needed
%                   when plotImage is set to 1. 
%
% Outputs:
%   - Im containing the MIP image. 

%if size(voxelSize,2) ~= 3 
if ~(sum(size(voxelSize) == [3 1]) || sum(size(voxelSize) == [1 3]))  
    resamp = 0;
    disp('No resampling is performed, since the parameter');
    disp('"voxelSize" does not contain a vector of length three.');
end


if resamp == 1
    vol = getIsotropicVoxels(vol, voxelSize(1), voxelSize(2), ...
                             voxelSize(3), 'cubic');   
end

im = maximumImageProjection(vol,dir);

if (plotImage == 1)
    if nargin == 6
        if figNumber > 0
            figure(round(figNumber));
        else
            figure;
        end
    else
        figure;
    end
    
    imagesc(im,[0 max(im(:))]); colormap gray; axis image; %colorbar
end