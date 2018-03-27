function [displacement, certaintyInc] = nonLinearRegistration(fixedVol,...
                                                              deformVol,...
                                                             kernelSize,...
                                                             tensorSigma)
                                                          

% Function that performs an iteration of the phase-based non-linear 
% registration. It calculates the incremental displacement field and the 
% incremental certainty for fitting the target volume (deformVol) to the 
% reference volume (fixedVol). The entire non-linear registration is given 
% by the 'morphon' function.
%
% Input:   fixedVol:    The reference volume that the target volume 
%                       (atlas) is registered to. 
%           
%          deformVol:   The volume that is to be registered to the
%                       reference volume.
%
%          kernelSize:  The size of the Gaussian kernel used for the 
%                       smoothing of A, b and the structure tensors.
%           
%          tensorSigma: The standard deviation of the Gaussian kernel used 
%                       for the smoothing of A, b and the structure tensor.
%       
% 
% Output:  displacement:    A 4D volume where the 3 first dimenstions
%                           containt the incremental displacement field and
%                           the 4th dimension is used for the directions y,
%                           x and z.
%
%          certaintyInc:    3D volume containing the incremental certainty.
%
%   See also MORPHON.

                                                          
% Load filters for the Morphon
load('quadratureFiltersForMorphonRegistration3D.mat')                               

targetSize  = size(fixedVol);

% Difference in local phase
dphi      = zeros([targetSize 6]);
% Certainty measurment
certainty = zeros([targetSize 6]);

% Set the tensor elements to zero
tensorElement11 = zeros(targetSize);
tensorElement12 = zeros(targetSize);
tensorElement13 = zeros(targetSize);
tensorElement22 = zeros(targetSize);
tensorElement23 = zeros(targetSize);
tensorElement33 = zeros(targetSize);





for k = 1:6
    % Filter responses of the volume that we deform
    qTarget = convn(fixedVol,f{k},'same');
    qDeform = convn(deformVol,f{k},'same');
    
    filterProd = qDeform.*conj(qTarget);
    
    % Difference in local phase
    % dphi(:,:,:,k) = angle(q_deform.*conj(q_target(:,:,:,k)));
    dphi(:,:,:,k) = atan2(imag(filterProd),real(filterProd));
    
    % Certainty
    certainty(:,:,:,k) = sqrt(abs(filterProd)).*(cos(dphi(:,:,:,k)/2).^2);
    
    % Calculate structure tensor components
    absTarget = abs(qTarget);
    
    tensorElement11 = tensorElement11 + absTarget*m11{k};   % y-dir
    tensorElement12 = tensorElement12 + absTarget*m12{k};
    tensorElement13 = tensorElement13 + absTarget*m13{k};
    tensorElement22 = tensorElement22 + absTarget*m22{k};   % x-dir
    tensorElement23 = tensorElement23 + absTarget*m23{k};
    tensorElement33 = tensorElement33 + absTarget*m33{k};   % z-dir
end
 
    

% Smooth the structure tensor components
tensorCert = sqrt(tensorElement11.^2 + 2*tensorElement12.^2 + ...
    2*tensorElement13.^2 + tensorElement22.^2   + ...
    2*tensorElement23.^2 + tensorElement33.^2);


tensorElement11 = normalizedAveraging(tensorElement11, kernelSize,...
                                      tensorSigma ,tensorCert);

tensorElement12 = normalizedAveraging(tensorElement12, kernelSize,...
                                      tensorSigma ,tensorCert);

tensorElement13 = normalizedAveraging(tensorElement13, kernelSize,...
                                      tensorSigma ,tensorCert);

tensorElement22 = normalizedAveraging(tensorElement22, kernelSize,...
                                      tensorSigma ,tensorCert);

tensorElement23 = normalizedAveraging(tensorElement23, kernelSize,...
                                      tensorSigma ,tensorCert);

tensorElement33 = normalizedAveraging(tensorElement33, kernelSize,...
                                      tensorSigma ,tensorCert);




% Normalize the tensor components
tensorNorm = sqrt(tensorElement11.^2 + 2*tensorElement12.^2  + ...
                  tensorElement22.^2 + 2*tensorElement13.^2  + ...
                  tensorElement33.^2 + 2*tensorElement23.^2) + eps;

tensor11 = tensorElement11./max(tensorNorm(:));
tensor12 = tensorElement12./max(tensorNorm(:));
tensor13 = tensorElement13./max(tensorNorm(:));
tensor22 = tensorElement22./max(tensorNorm(:));
tensor23 = tensorElement23./max(tensorNorm(:));
tensor33 = tensorElement33./max(tensorNorm(:));
    

% Calculation of the incremental displacement field and certainty 
[displacement, certaintyInc] = calculateDisplacement(tensor11, tensor12,...
                                       tensor13, tensor22, tensor23, ...
                                       tensor33, certainty, ...
                                       filterDirection, dphi);  


% Scale independent amplification factor for the displacement field
amplification = 2;

displacement(:,:,:,2) = amplification*displacement(:,:,:,2);
displacement(:,:,:,1) = amplification*displacement(:,:,:,1); 
displacement(:,:,:,3) = amplification*displacement(:,:,:,3); 


 

end