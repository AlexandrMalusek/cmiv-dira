function [finalDeformedVolume, dispAcc] = morphon(fixedVolume,   ...
                                                  deformedVolume,...
                                                  scales,        ...
                                                  iterations)
                                        
% Function for performing non-linear registration in 3D using the Morphon 
% algorithm. The Morphon uses the local phase of quadrature filters to fit 
% lines and edges in the target volume ('deformedVolume') to lines and 
% edges in the reference volume ('fixedVolume') respectively.
%
% 
% Inputs:   - fixedVolume:      The volume that the deformedVolume is 
%                               fitted to by the Morphon algorithm.
%
%           - deformedVolume:   The volume that should be deformed in order
%                               to fit the fixed volume (reference volume).
%                               This volume can be an atlas or another data
%                               set.
%
%           - scales:           Vector containing the factor of the highest
%                               downsampling (coarsest scale) and the
%                               lowest downsampling (finest scale). The
%                               vector [5, 3] would mean that the
%                               registration would be done when the volumes
%                               have been downsampled by a factor 2^(-5), 
%                               2^(-4) and 2^(-3). 
%
%           - iterations:       Vector containing the number of iterations
%                               on each scale. The first entry is the
%                               number of iterations on the original scale
%                               (downsampling 2^(-0)), the next entry is  
%                               for the number of iterations on the second 
%                               finest scale (downsampling 2^(-1)), ...
%
%
% Outputs:  
%       - finalDeformedVolume:  The volume defomed by the Morphon
%                               algorithm.
%
%       - dispAcc:              The accumulated displacement field. The 4th
%                               dimension is used for the y, x and z
%                               dimension.
%
%
%
%
% The implementation of the Morphon algorithm is based on the following
% articles:
%
%   1. "MORPHONS : SEGMENTATION USING ELASTIC CANVAS AND PAINT ON PRIORS"
%
%   2. "Phase-Based Non-Rigid 3D Image Registration: From Minutes to 
%       Seconds Using CUDA"


% %%%%% ------ Parameters ------ %%%%% %
% Regularization of accumulative displacement field
filterSizeAcc = [13, 13,  25, 25, 13,   9,   9]; 
sigmaAcc      = [3.0,  3.0,   2.0, 1.8, 1.1,  0.7,   2]; 

% Smoothing of A, b and the structure tensors when calculating the 
% incremental displacement field and certainty
filterSizeInc = [13 13 11 7 7 7 7 7];
sigmaInc      = [2.5 2.5 2 1 1 1 1 1];



% %%%%% ------ Initialization ------ %%%%% %

fixedSize = size(fixedVolume);

% Accumulated displacement
dispAcc = zeros([ceil(fixedSize*2^(-scales(1))) 3]);

% Accumulated certainty 
certAcc = zeros([ceil(fixedSize*2^(-scales(1))) 3]);

% Rescale the volume that plan to deform
if (scales(1) > 0)
    smallDeformed  = rescale3(deformedVolume,                ...
                             ceil(fixedSize*2^(-scales(1))), ...
                             'linear');
                         
    downsampledDeformed = smallDeformed;
else
    smallDeformed  = deformedVolume;
    downsampledDeformed = smallDeformed;
end



% %%%%% ------ Main loop ------ %%%%% %

for rho = scales(1):-1:scales(2) % rho is the scale parameter 
    
    % Calculate the size of the new downsampled image
    smallSize = ceil(fixedSize*2^( -rho));
    nextSize  = ceil(fixedSize*2^(1 - rho));

    [y,x,z]   = ndgrid(1:smallSize(1),1:smallSize(2), 1:smallSize(3));

    % Resizing of the fixed volume (reference volume) 
    if(rho > 0)
        smallFixed     = rescale3(fixedVolume,smallSize,'linear'); 
    else 
        smallFixed     = fixedVolume;
    end
       
    for iter = 1:iterations(rho+1)
        
       
        kernelSize  = filterSizeInc(rho+1);
        tensorSigma = sigmaInc(rho+1);


        [dispInc, certInc] = nonLinearRegistration(smallFixed,    ...
                                                   smallDeformed, ...
                                                   kernelSize,    ...
                                                   tensorSigma);
                                               
        certInc3 = repmat(certInc,1,1,1,3);
        
        % Update the accumulated displacement field
        dispAcc = (certAcc.*dispAcc + certInc3.*(dispAcc + dispInc))...
               ./ (certAcc + certInc3 + eps);
        
        % Update the accumulated certainty
        certAcc = (certAcc.^2 + (2^(-rho)* certInc3).^2)./ ...                                                                               %% CANVAS AND PAINT ON PRIORS
                  (certAcc    +  2^(-rho)* certInc3 + eps);

        for m = 1:3
             dispAcc(:,:,:,m) = normalizedAveraging(dispAcc(:,:,:,m),...
                                               filterSizeAcc(rho+1), ...
                                               sigmaAcc(rho+1), ...
                                               certAcc(:,:,:,m));
        end

        % Upsampling of the accumulated certainty and displacement field to
        % the next scale. 
        if(rho > scales(2) && iter == iterations(rho+1)) 

            [dispAcc, certAcc] = changeOfScale(dispAcc, certAcc, ...
                                               smallSize, nextSize);
                                              
            [y,x,z] = ndgrid(1:nextSize(1),1:nextSize(2), 1:nextSize(3));
                                  
                                              
            % downsampling of the original atlas to the size 'nextSize'.
            downsampledDeformed  = rescale3(deformedVolume, ...
                                            nextSize,       ...
                                            'linear');                                              
                                        
        elseif(rho==scales(2) && iter == iterations(rho+1))
            [dispAcc, certAcc] = changeOfScale(dispAcc, certAcc,    ...
                                               smallSize, fixedSize);
                                           
            [y,x,z] = ndgrid(1:fixedSize(1),1:fixedSize(2),1:fixedSize(3));
            
            smallFixed = fixedVolume;
            downsampledDeformed  = deformedVolume;
        elseif(rho == 0)
                downsampledDeformed  = deformedVolume;               
        end
 
       % Apply the displacement to the deformed volume (target volume)                     
       smallDeformed = ba_interp3(x,y,z,double(downsampledDeformed),    ...
                                        x - double(dispAcc(:,:,:,2)),   ...
                                        y - double(dispAcc(:,:,:,1)),   ...
                                        z - double(dispAcc(:,:,:,3)));                                            
       
       smallDeformed(isnan(smallDeformed)) = 0;

    end
end

finalDeformedVolume = smallDeformed;

end