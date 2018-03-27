function [transformedAtlas, ...
          displacement, pScales] = linearRegistration(volume, atlas,  ...
                                          scales, iterations, weight, ...          
                                          initialTranslationXY)
% Function that fits an atlas to a volume using phase-based affine 
% registration. The linear registration is performed over different scales,
% first fitting the atlas on a coarse scale to the volume and then refining 
% the fitting on the finer scales. 
%
% It is good if the atlas covers a larger section of the body than the 
% volume it is registered to in order to avoid getting regions with missing
% information. This function allows to register an altas that contains more
% image slices than the volume to the volume. This is done by calculating   
% the transformation for the slices in the center of the atlas (same amount 
% of slices as the volume) and thereafter apply the found transformation to 
% the entire atlas. This is done in each iteration of the algorithm, which 
% allows the linear registration to move parts of the image content from 
% outside the central region (central slizes) inside the central region or
% the other way arround. 
%
%
% Input:  - volume:     3D volume.
%
%         - atlas:      Atlas to fit to the 3D volume.
%
%         - scales:     A vector containing the coarsest scale and the 
%                       finest scale. Both are described as integers, where
%                       scale 0 means the original scale (no downsampling).
%
%         - iterations: Vector with the number of iterations on each of the
%                       scales. The first entry is the number of iterations 
%                       on the original scale, the next entry is for when 
%                       both volumes are downsampled with a factor 2 and so
%                       on.
%
%          - weight:    Weight function that is multiplied to the certainty
%                       measurement c, when calculation the parameter 
%                       vector. 
%
%          - initialTranslationXY:
%                       A vector containing a translation in x and y 
%                       direction that is used as initial translation for
%                       the parameter vector. 
%
%
% Output:  - transformedAtlas:  The transformed atlas.
%                       
%          - displacement:      The final displacement field. It is
%                               estimated for the original size of the
%                               atlas volume. 
%
%          -  pScales:          The parameter vector.


%% Initialize varaibles
volSize   = size(volume);
atlasSize = size(atlas);

scaleStart = scales(1);
scalesEnd  = scales(2);  
pScales    = zeros(12,1);

% Scaling the atlas to the coarsest scale
smallAtlas = rescale3(atlas,ceil(2^(-scaleStart)*atlasSize),'cubic');            

%% Initial placement of the atlas 
% Perform the initial translation in X and Y direction
pScales(1:2) = initialTranslationXY*2^(-scaleStart);

% Place the atlas in the middle of the volume
pScales(3) = max((atlasSize(3)-volSize(3))/2,0)*2^( -scaleStart);

smallSizeA = ceil(atlasSize*2^( -scaleStart));

motionVectorX = ones(smallSizeA)*pScales(1);
motionVectorY = ones(smallSizeA)*pScales(2);
motionVectorZ = ones(smallSizeA)*pScales(3);

[y,x,z] = ndgrid(-floor((smallSizeA(1)-1)/2):ceil((smallSizeA(1)-1)/2), ...       
                 -floor((smallSizeA(2)-1)/2):ceil((smallSizeA(2)-1)/2), ...
                 -floor((smallSizeA(3)-1)/2):ceil((smallSizeA(3)-1)/2));

% Get the atlas moved by the initial translation 
smallAtlas = ba_interp3(x,y,z,double(smallAtlas), ...
                         x + motionVectorX, ...
                         y + motionVectorY, ...
                         z + motionVectorZ, 'linear');
smallAtlas(isnan(smallAtlas)) = 0;

%% Main loop of the affine registration

for rho = scaleStart:-1:scalesEnd 
    
    % Calculate the size of the new downsampled image
    smallSizeV      = ceil(volSize  *2^(-rho));
    nextSizeV       = ceil(volSize  *2^(1 - rho));
    nextSizeAtlas   = ceil(atlasSize*2^(1 - rho));
    
    % Resize our reference volume
    smallVolume   = rescale3(volume, smallSizeV,'cubic');                          
    smallWeight   = rescale3(weight, smallSizeV,'cubic');                          
       
    [~ , ~, p] = getAffineTransformation(smallVolume, smallAtlas, ...
                                         iterations(rho+1),smallWeight);
                                     
    clear small_fixed small_weight; 
    
    % update the the displacement field
    pScales = pScales + p;

    if(rho > scalesEnd) % If there are finer scales
       lengthDiffInZ = nextSizeAtlas(3) - nextSizeV(3);
       
       [y,x,z] = ndgrid(-floor((nextSizeV(1)-1)/2):ceil((nextSizeV(1)-1)/2),...       
                        -floor((nextSizeV(2)-1)/2):ceil((nextSizeV(2)-1)/2),...
                        -floor((nextSizeV(3)-1)/2):ceil((nextSizeV(3)-1)/2) ...
                        + lengthDiffInZ);
                 
       % Rescale the translation      
       scale = nextSizeV./smallSizeV;
       pScales(1:3) =  scale'.*pScales(1:3);
       motionVectorX = ones(nextSizeAtlas)*pScales(1) + x*pScales(4) +...
                                         y*pScales(5) + z*pScales(6);      
       motionVectorY = ones(nextSizeAtlas)*pScales(2) + x*pScales(7) +...
                                         y*pScales(8) + z*pScales(9);
       motionVectorZ = ones(nextSizeAtlas)*pScales(3) + x*pScales(10)+...
                                         y*pScales(11)+ z*pScales(12);   

 
       % Downsampling of the atlas to the next size. 
       downsampledAtlas = rescale3(atlas,nextSizeAtlas,'cubic');  
       
       % Apply the transformation to the downsampled atlas.
       smallAtlas = ba_interp3(x,y,z,double(downsampledAtlas), ...
                                        x + motionVectorX,   ...
                                        y + motionVectorY,   ...
                                        z + motionVectorZ);                                       
                                    
                                
       smallAtlas(isnan(smallAtlas)) = 0;
       smallAtlas = single(smallAtlas);
       clear downsampledAtlas motionVectorX motionVectorY motionVectorZ;  
        

    else % If we are at the finest scale
       lengthDiffInZ = atlasSize(3) - volSize(3);
       [y,x,z] = ndgrid(-floor((volSize(1)-1)/2):ceil((volSize(1)-1)/2),...      
                        -floor((volSize(2)-1)/2):ceil((volSize(2)-1)/2),...
                        -floor((volSize(3)-1)/2):ceil((volSize(3)-1)/2) ...
                        + lengthDiffInZ);

       % Rescale the translation
       scale = volSize./smallSizeV;
       pScales(1:3) =  scale'.*pScales(1:3);   
       
       % Calculate motion vector
       motionVectorX = ones(atlasSize)*pScales(1)  + x*pScales(4)  + ...
                                     y*pScales(5)  + z*pScales(6);            
       motionVectorY = ones(atlasSize)*pScales(2)  + x*pScales(7)  + ...
                                     y*pScales(8)  + z*pScales(9);
       motionVectorZ = ones(atlasSize)*pScales(3)  + x*pScales(10) + ...
                                     y*pScales(11) + z*pScales(12);   
       
       % Apply transformation on the atlas (original scale) 
       transformedAtlas = ba_interp3(x,y,z,double(atlas), ...
                                    x + motionVectorX,   ...
                                    y + motionVectorY,   ...
                                    z + motionVectorZ,'linear');                          
       transformedAtlas(isnan(transformedAtlas)) = 0;
       
    end
end

% Save the final displacement field.
displacement = zeros([atlasSize 3]);
displacement(:,:,:,1) = motionVectorY;
displacement(:,:,:,2) = motionVectorX;
displacement(:,:,:,3) = motionVectorZ;

end