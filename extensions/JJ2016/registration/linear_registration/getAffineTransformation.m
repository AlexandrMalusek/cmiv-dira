function [registeredAtlas, ...
          displacement,... 
          pTotal]   =  getAffineTransformation(volume, atlas, ...
                                               iterations,certaintyWeight) 
                                           
% Function that fits an atlas to an volume by an affine transformation. The
% registration is based on the local phase which fits edges to edges and 
% lines to lines. 
%
% The function allows that the used atlas has more elements in z direction 
% than the volume the atlas is registered to. The difference in voxel array
% dimensionallity is handled by using the k first elements in z direction 
% of the atlas for the estimation of the parameter vector, where k is the
% number of elements of the volume in z direction. A displacement field is 
% thereafter calculated for the entire atlas volume. This is repeated in
% each iteration which allows that the entire content of the atlas can be
% used in the registration.
%
% Input:    - volume            The 3D volume to which the atlas shall be 
%                               fitted.
%
%           - atlas             The 3D atlas.
%
%           - iterations        The number of iterations used to fit the 
%                               atlas to the volume.
%
%           - certaintyWeight   Volume with the same size as the volume. It
%                               can be used to reduce the certainty of the 
%                               voxels near the boarders of the volumes.
%
%
%
% Output:   - registeredAtlas   The registered atlas volume.
%
%           - displacement      The calculated displacement field
%                               describing the affine transformation. 
%
%           - pTotal            The parameter vector describing the affine
%                               transformation.
%
%   See also linearRegistration.


% Get the filter responses of the volume
load('filters_for_parametric_registration.mat');

qVol1 = convn(volume,f1_parametric_registration,'same'); 
qVol2 = convn(volume,f2_parametric_registration,'same');
qVol3 = convn(volume,f3_parametric_registration,'same');


% ----- Initialize variables ----- %
[sizeY,  sizeX,  sizeZ ] = size(volume);   
[atlasY, atlasX, atlasZ] = size(atlas);     

% Create a meshgrid                      
[y,x,z] = ndgrid(-floor((sizeY-1)/2) : ceil((sizeY-1)/2),...
                 -floor((sizeX-1)/2) : ceil((sizeX-1)/2),...
                 -floor((sizeZ-1)/2) : ceil((sizeZ-1)/2));            
             
% Create a meshgrid                             
lengthDiffZ = atlasZ - sizeZ;

[ay,ax,az] = ndgrid(-floor((sizeY-1)/2) : ceil((sizeY-1)/2),...     
                    -floor((sizeX-1)/2) : ceil((sizeX-1)/2),...
                    -floor((sizeZ-1)/2) : ceil((sizeZ-1)/2) + lengthDiffZ);                        

% Total Parameter vector
pTotal    = zeros(12,1);

% Motion vectors
motionVectorY = zeros(atlasY, atlasX, atlasZ);
motionVectorX = zeros(atlasY, atlasX, atlasZ);     
motionVectorZ = zeros(atlasY, atlasX, atlasZ);

% Gradients
gradX = zeros(sizeY, sizeX, sizeZ);
gradY = zeros(sizeY, sizeX, sizeZ);       
gradZ = zeros(sizeY, sizeX, sizeZ);


% ----- Main Loop ----- %

for iter = 1 : iterations
    atlasVolume = ba_interp3(ax,ay,az,double(atlas),ax + motionVectorX, ...
                                                    ay + motionVectorY, ...
                                                    az + motionVectorZ);                     
                                         
    atlasVolume(isnan(atlasVolume)) = 0;
    
    % crop atlas here
    croppedAtlas = atlasVolume(1:sizeY,1:sizeX, 1:sizeZ);
    
    % Filter responses of 
    qRef1 = convn(croppedAtlas,f1_parametric_registration,'same'); 
    qRef2 = convn(croppedAtlas,f2_parametric_registration,'same');
    qRef3 = convn(croppedAtlas,f3_parametric_registration,'same');
    
    % Phase
    dphi1 = angle(qVol1.*conj(qRef1));
    dphi2 = angle(qVol2.*conj(qRef2));
    dphi3 = angle(qVol3.*conj(qRef3));
    
    % Certainty
    c1 = sqrt(abs(qRef1.*qVol1)).*((cos(dphi1/2)).^2);
    c2 = sqrt(abs(qRef2.*qVol2)).*((cos(dphi2/2)).^2);
    c3 = sqrt(abs(qRef3.*qVol3)).*((cos(dphi3/2)).^2);    

    % Decreases the certainty near the edges of the volume.
    c1 = certaintyWeight .* c1; 
    c2 = certaintyWeight .* c2;
    c3 = certaintyWeight .* c3;
    
    % Gradient
    gradX(:,2:end-1,:) = angle(qVol1(:,3:end,:).*conj(qVol1(:,2:end-1,:))...
                    +  qVol1(:,2:end-1,:).*conj(qVol1(:,1:end-2,:))  ...
                    +  qRef1(:,3:end,:)  .*conj(qRef1(:,2:end-1,:))  ...
                    +  qRef1(:,2:end-1,:).*conj(qRef1(:,1:end-2,:)));
    
    gradY(2:end-1,:,:) = angle(qVol2(3:end,:,:).*conj(qVol2(2:end-1,:,:))...
                    +  qVol2(2:end-1,:,:) .*conj(qVol2(1:end-2,:,:)) ...
                    +  qRef2(3:end,  :,:) .*conj(qRef2(2:end-1,:,:)) ...
                    +  qRef2(2:end-1,:,:) .*conj(qRef2(1:end-2,:,:)));
                   
    gradZ(:,:,2:end-1) = angle(qVol3(:,:,3:end).*conj(qVol3(:,:,2:end-1))...
                    +  qVol3(:,:,2:end-1) .*conj(qVol3(:,:,1:end-2)) ...
                    +  qRef3(:,:,3:end  ) .*conj(qRef3(:,:,2:end-1)) ...
                    +  qRef3(:,:,2:end-1) .*conj(qRef3(:,:,1:end-2)));
     
    
    % Calculate h
    h = [sum(c1(:).^2.*dphi1(:).*gradX(:)); ...
         sum(c2(:).^2.*dphi2(:).*gradY(:)); ...
         sum(c3(:).^2.*dphi3(:).*gradZ(:)); ...
         sum(x(:).*c1(:).^2.*dphi1(:).*gradX(:));...
         sum(y(:).*c1(:).^2.*dphi1(:).*gradX(:));...
         sum(z(:).*c1(:).^2.*dphi1(:).*gradX(:));...
         sum(x(:).*c2(:).^2.*dphi2(:).*gradY(:));...
         sum(y(:).*c2(:).^2.*dphi2(:).*gradY(:));...
         sum(z(:).*c2(:).^2.*dphi2(:).*gradY(:));...
         sum(x(:).*c3(:).^2.*dphi3(:).*gradZ(:));...
         sum(y(:).*c3(:).^2.*dphi3(:).*gradZ(:));...
         sum(z(:).*c3(:).^2.*dphi3(:).*gradZ(:))];

    % Calculate A 
    A = zeros(12,12);

    A(1,1)  = sum(c1(:).^2 .*gradX(:).^2);
    A(2,2)  = sum(c2(:).^2 .*gradY(:).^2);
    A(3,3)  = sum(c3(:).^2 .*gradZ(:).^2);
    
    A(4,4)  = sum(c1(:).^2 .*x(:).^2 .*gradX(:).^2);
    A(5,5)  = sum(c1(:).^2 .*y(:).^2 .*gradX(:).^2);
    A(6,6)  = sum(c1(:).^2 .*z(:).^2 .*gradX(:).^2);    

    A(7,7)  = sum(c2(:).^2 .*x(:).^2 .*gradY(:).^2);
    A(8,8)  = sum(c2(:).^2 .*y(:).^2 .*gradY(:).^2);
    A(9,9)  = sum(c2(:).^2 .*z(:).^2 .*gradY(:).^2);
    
    A(10,10) = sum(c3(:).^2 .*x(:).^2 .*gradZ(:).^2);
    A(11,11) = sum(c3(:).^2 .*y(:).^2 .*gradZ(:).^2);
    A(12,12) = sum(c3(:).^2 .*z(:).^2 .*gradZ(:).^2);    


    
    A(1,4) = sum(c1(:).^2 .* x(:).*gradX(:).^2);
    A(1,5) = sum(c1(:).^2 .* y(:).*gradX(:).^2);
    A(1,6) = sum(c1(:).^2 .* z(:).*gradX(:).^2);
    
    A(2,7) = sum(c2(:).^2 .* x(:) .*gradY(:).^2);
    A(2,8) = sum(c2(:).^2 .* y(:) .*gradY(:).^2);
    A(2,9) = sum(c2(:).^2 .* z(:) .*gradY(:).^2);
    
    A(3,10) = sum(c3(:).^2 .* x(:) .*gradZ(:).^2);
    A(3,11) = sum(c3(:).^2 .* y(:) .*gradZ(:).^2);
    A(3,12) = sum(c3(:).^2 .* z(:) .*gradZ(:).^2);

    
    A(4,5) = sum(c1(:).^2 .* x(:) .*y(:) .*gradX(:).^2);
    A(4,6) = sum(c1(:).^2 .* x(:) .*z(:) .*gradX(:).^2);
    A(5,6) = sum(c1(:).^2 .* y(:) .*z(:) .*gradX(:).^2);    

    A(7,8) = sum(c2(:).^2 .* x(:) .*y(:) .*gradY(:).^2);
    A(7,9) = sum(c2(:).^2 .* x(:) .*z(:) .*gradY(:).^2);
    A(8,9) = sum(c2(:).^2 .* y(:) .*z(:) .*gradY(:).^2);    
    
    A(10,11) = sum(c3(:).^2 .* x(:) .*y(:) .*gradZ(:).^2);
    A(10,12) = sum(c3(:).^2 .* x(:) .*z(:) .*gradZ(:).^2);
    A(11,12) = sum(c3(:).^2 .* y(:) .*z(:) .*gradZ(:).^2);      
     
    Atot=A+triu(A,1)';
    
    
    % Get p
    p       = Atot\h;
    pTotal  = pTotal+p;

    pTotal  = double(pTotal);
    
    onesAtlas = ones(atlasY, atlasX, atlasZ);
    % Motion Vector
    motionVectorY = onesAtlas*pTotal(2)  + ax*pTotal(7)  + ...
                           ay*pTotal(8)  + az*pTotal(9);
                       
    motionVectorX = onesAtlas*pTotal(1)  + ax*pTotal(4)  + ...
                           ay*pTotal(5)  + az*pTotal(6);
                       
    motionVectorZ = onesAtlas*pTotal(3)  + ax*pTotal(10) + ...
                           ay*pTotal(11) + az*pTotal(12);
   
                       
    motionVectorY = double(motionVectorY);
    motionVectorX = double(motionVectorX);    
    motionVectorZ = double(motionVectorZ);   
end

% Output parameter
registeredAtlas = ba_interp3(ax,ay,az,double(atlas),...
                                  ax + motionVectorX, ...
                                  ay + motionVectorY, ...
                                  az + motionVectorZ);                     
                                         
registeredAtlas(isnan(registeredAtlas)) = 0;

displacement = zeros(atlasY, atlasX, atlasZ,3);
displacement(:,:,:,1) = motionVectorX;
displacement(:,:,:,2) = motionVectorY;
displacement(:,:,:,3) = motionVectorZ;




end