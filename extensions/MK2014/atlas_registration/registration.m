
function [vx,vy] = registration(image, atlas, num_scales,iter)
% Registration function. This function performs an affine transformation in
% order to match the atlas (second input) on the image (first input). The
% registration starts on a coarse scale and is thereafter updated by
% using finer scales.  
%
% Input: Patient image:    The gray scale image that the atlas is matched 
%                          to. Needs to be double or single.
%
%        Atlas image:      Image that will be matched on the patient image.
%                          Needs to be double or single.
%
%        Number of scales: The number of scales used for the registration. 
%                          The downsampling is done by a factor of two.
%
%        Iterations:       A vector with the number of iterations on each 
%                          scale. The first element of the vector is 
%                          corresponding to the number of iterations on the
%                          finest scale (no downsampling) and element k 
%                          corresponds to the number of iterations on scale 
%                          k-1.

% --------------------------------------------------------------------------
% Create Scale-Pyramid------------------------------------------------------

pyramidImage = cell(1,num_scales);
pyramidAtlas = cell(1,num_scales);
scale        = cell(1,num_scales);

size_y = size(image,1);

for k = 0:num_scales
    scale{k+1} = min(size_y*2^(-k) +1,size_y)/size_y;
    
    pyramidImage{k+1} = double(imresize(image, scale{k+1}));
    pyramidAtlas{k+1} = double(imresize(atlas, scale{k+1}));
    
end

% --------------------------------------------------------------------------
% Registration from coarse to fine scale------------------------------------

smallAtlasMoved = pyramidAtlas{num_scales+1};
ptot = zeros(6,1);

for k = (num_scales+1):-1:2

    p = registerAtlas(pyramidImage{k},smallAtlasMoved,iter(k));

    %create a grid for the small image
    [sy_old, sx_old] = size(pyramidAtlas{k});
    [sy, sx] = size(pyramidAtlas{k-1});
    [x,y] = meshgrid(-(sx-1)/2:(sx-1)/2,-(sy-1)/2:(sy-1)/2);
    
    if (k>2)
        % Rescaling the translation
        scale_factor = [sy/sy_old, sx/sx_old];
        p(1:2) = scale_factor'.*p(1:2);
        ptot = ptot + p;
        
        vx=1*ptot(1)+x*ptot(3)+y*ptot(4);
        vy=1*ptot(2)+x*ptot(5)+y*ptot(6);
        
        % Upsamling of Motion Field
    
        smallAtlasMoved = interp2(x,y,double(pyramidAtlas{k-1}),x+vx,y+vy);
        smallAtlasMoved(isnan(smallAtlasMoved)) = 0;
        
    else
        % In the last Scale we don't need to rescale the translation
        ptot = ptot + p;
        
        vx=1*ptot(1)+x*ptot(3)+y*ptot(4);
        vy=1*ptot(2)+x*ptot(5)+y*ptot(6);
    end 
    
end

end
