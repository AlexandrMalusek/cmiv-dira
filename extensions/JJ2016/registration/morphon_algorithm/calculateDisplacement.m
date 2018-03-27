function [displacement, certaintyInc] = calculateDisplacement(t11,t12,...
                                                    t13,t22,t23,t33,...
                                                    certain,filterDir,dphi)
% Function that calculates the incremental displacement field and the 
% incremental certainty. This function is a part of the implementation of 
% the Morphon algorithm, see 'morphon' function. 
%
% Inputs:   - The elements of the tensor for the voxels given as the  
%             following 6 volumes: t11, t12, t13, t22, t23, t33.
%
%           - A volume containg the certainty for the 6 filter directions.
%
%           - A cell array containing the filter directions of the 6
%             quadrature filters.
%
%           - The difference in local phase.
%
%
%
% Output:   - The incremental displacement field in x, y and z direction. 
%             The 4th dimension is used for the direction of the 
%             displacement field. They are ordered y, x and z. 
%
%           - The incremental certainty.
% 
%
% Example: 
% 
%  [displacement, certaintyInc] = calculateDisplacement(t11,t12,...
%                                                   t13,t22,t23,t33,...
%                                                   certain,filterDir,dphi)
%
%   See also Morphon.

sigma = 1.5;                                           

tensorSize = size(t11);                                             
% Initialize the elements of the A matrix and the b vector                                             
A11 = zeros(tensorSize);
A12 = zeros(tensorSize);
A13 = zeros(tensorSize);
A22 = zeros(tensorSize);
A23 = zeros(tensorSize);
A33 = zeros(tensorSize);

b1 = zeros(tensorSize);
b2 = zeros(tensorSize);
b3 = zeros(tensorSize);

% Calculate T'* T
tensor11 = t11.^2 + t12.^2 + t13.^2;
tensor12 = t11.*t12 + t12.*t22 + t13.*t23;
tensor13 = t11.*t13 + t12.*t23 + t13.*t33;
tensor22 = t12.^2 +t22.^2 + t23.^2;
tensor23 = t12.*t13+ t22.*t23 + t23.*t33;
tensor33 = t13.^2 + t23.^2 + t33.^2;

%
for k = 1:6
    % Calculate A 
    A11 = A11 + certain(:,:,:,k).*tensor11; 
    A12 = A12 + certain(:,:,:,k).*tensor12; 
    A13 = A13 + certain(:,:,:,k).*tensor13; 
    A22 = A22 + certain(:,:,:,k).*tensor22; 
    A23 = A23 + certain(:,:,:,k).*tensor23; 
    A33 = A33 + certain(:,:,:,k).*tensor33; 
    
    % Calculate b
    b1 = b1+certain(:,:,:,k).*dphi(:,:,:,k).*(filterDir{k}(1)*tensor11 ...
                                            + filterDir{k}(2)*tensor12 ...
                                            + filterDir{k}(3)*tensor13); 
                                      
    b2 = b2+certain(:,:,:,k).*dphi(:,:,:,k).*(filterDir{k}(1)*tensor12 ...
                                           +  filterDir{k}(2)*tensor22 ...
                                           +  filterDir{k}(3)*tensor23);     
                                      
    b3 = b3+certain(:,:,:,k).*dphi(:,:,:,k).*(filterDir{k}(1)*tensor13 ...
                                           +  filterDir{k}(2)*tensor23 ...
                                           +  filterDir{k}(3)*tensor33);                            
end

% Smooth the components of A and b
A11 = gaussSmoothing(A11,sigma);
A12 = gaussSmoothing(A12,sigma);
A13 = gaussSmoothing(A13,sigma);
A22 = gaussSmoothing(A22,sigma);
A23 = gaussSmoothing(A23,sigma);
A33 = gaussSmoothing(A33,sigma);

b1 = gaussSmoothing(b1,sigma);
b2 = gaussSmoothing(b2,sigma);
b3 = gaussSmoothing(b3,sigma);

% Calculate certainty 
certaintyInc = A11 + A22 + A33;   

norm = 1./(A11.*(A22.*A33 - A23.^2)   + A12.*(A23.*A13 - A12.*A33) ...
             + A13.*(A12.*A23 - A22.*A13) + eps);  
    
displacement(:,:,:,1) = norm.*((A22.*A33-A23.*A23).*b1 ...
                           + (A13.*A23-A12.*A33).*b2 ...  
                           + (A12.*A23-A13.*A22).*b3);   
    
displacement(:,:,:,2)  = norm.*((A23.*A13-A12.*A33).*b1 ...
                           + (A11.*A33-A13.*A13).*b2 ...
                           + (A13.*A12-A11.*A23).*b3);
    
                       
displacement(:,:,:,3) = norm.*((A12.*A23-A13.*A22).*b1 ...
                           + (A12.*A13-A11.*A23).*b2 ...
                           + (A11.*A22-A12.*A12).*b3);

displacement(isnan(displacement)) = 0; 


end

