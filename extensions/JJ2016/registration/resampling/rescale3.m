function resampledVolume = rescale3(volume, newSize, method)
% Function that resamples a 3D volume to the size given by newSize. 
%
% Input:    - volume:   Volume that is to be resampled.
%
%           - newSize:  The new size of the volume defined as a vector
%                       containing the new size in y, x and z direction.
%
%           - method:   The type of resampling used for the interpolation. 
%                       method can be set to 'nearest', 'linear' and 
%                       'cubic'. 
%
%
% Output:   - resampledVolume: The resampled volume.

% downsamp = faktor of downsampling, i.e. faktor 2
downsamp = size(volume)./newSize;

[sy,sx,sz]  = size(volume);

isOddSize   = ceil(mod(downsamp,2));
kernelSize  = 1 + ceil(downsamp) + isOddSize;
sigmaScaled = kernelSize/4;

if (downsamp(1) > 1)
    gaussKernelY = fspecial('Gaussian',[1 kernelSize(1)],sigmaScaled(1));
    volume = convn(volume, gaussKernelY', 'same');
end

if (downsamp(2) > 1)
    gaussKernelX = fspecial('Gaussian',[1 kernelSize(2)],sigmaScaled(2));
    volume = convn(volume, gaussKernelX , 'same');
end

if (downsamp(3) > 1)
    gaussKernelZ = fspecial('Gaussian',[1 kernelSize(3)],sigmaScaled(3));
    volume = convn(volume, permute(gaussKernelZ, [3 1 2]), 'same');
end

[y, x, z] = ndgrid(linspace(1,sy,newSize(1)),...        % y - size
                   linspace(1,sx,newSize(2)),...        % x - size
                   linspace(1,sz,newSize(3)));          % z - size

resampledVolume = ba_interp3(double(volume),x, y, z,method);
resampledVolume = single(resampledVolume);

% Could be used if ba_interp3 is not available. 
%resampledVolume = interp3(volume,x, y, z,Method); 
end