function [A11,A12,A13,A22,A23,A33,h1,h2,h3] = caclulateAAndB(t11,t12,t13,t22,t23,t33,...
                                                 certain,filterDir,d_k)
sigma = 1.5;                                           

tensorSize = size(t11);                                             
% Initialize the elements of the A matrix and the h vector                                             
A11 = zeros(tensorSize);
A12 = zeros(tensorSize);
A13 = zeros(tensorSize);
A22 = zeros(tensorSize);
A23 = zeros(tensorSize);
A33 = zeros(tensorSize);

h1 = zeros(tensorSize);
h2 = zeros(tensorSize);
h3 = zeros(tensorSize);

% Calculate T'* T
tensor11 = t11.^2 + t12.^2 + t13.^2;
tensor12 = t11.*t12 + t12.*t22 + t13.*t23;
tensor13 = t11.*t13 + t12.*t23 + t13.*t33;
tensor22 = t12.^2 +t22.^2 + t23.^2;
tensor23 = t12.*t13+ t22.*t23 + t23.*t33;
tensor33 = t13.^2 + t23.^2 + t33.^2;

%%
for k = 1:6
    % Calculate A 
    A11 = A11 + certain(:,:,:,k).*tensor11; 
    A12 = A12 + certain(:,:,:,k).*tensor12; 
    A13 = A13 + certain(:,:,:,k).*tensor13; 
    A22 = A22 + certain(:,:,:,k).*tensor22; 
    A23 = A23 + certain(:,:,:,k).*tensor23; 
    A33 = A33 + certain(:,:,:,k).*tensor33; 
    
    % Calculate h
    h1 = h1 + certain(:,:,:,k) .* d_k(:,:,:,k) .* (filterDir{k}(1)*tensor11 ...
                                                 + filterDir{k}(2)*tensor12 ...
                                                 + filterDir{k}(3)*tensor13); 
                                      
    h2 = h2 + certain(:,:,:,k) .* d_k(:,:,:,k) .* (filterDir{k}(1)*tensor12 ...
                                                +  filterDir{k}(2)*tensor22 ...
                                                +  filterDir{k}(3)*tensor23);     
                                      
    h3 = h3 + certain(:,:,:,k) .* d_k(:,:,:,k) .* (filterDir{k}(1)*tensor13 ...
                                                +  filterDir{k}(2)*tensor23 ...
                                                +  filterDir{k}(3)*tensor33);                            

end


A11 = gaussSmoothing(A11,sigma);
A12 = gaussSmoothing(A12,sigma);
A13 = gaussSmoothing(A13,sigma);
A22 = gaussSmoothing(A22,sigma);
A23 = gaussSmoothing(A23,sigma);
A33 = gaussSmoothing(A33,sigma);

h1 = gaussSmoothing(h1,sigma);
h2 = gaussSmoothing(h2,sigma);
h3 = gaussSmoothing(h3,sigma);


end