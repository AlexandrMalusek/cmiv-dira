function newIm = VisualizeSegmentationGray(im,mask,tissue)
% Redefine the values of an image inside a region defined by the mask. 
% The redefined values are the same as for the atlas and that result in the
% same gray levels in the image.  
%
% Input:    - image: Segmented gray scale image
%
%           - segmentation result: binary mask showing the segmentation 
%             result
%
%           - name of the tissue: This can be set to 'bone', 'prostate', 
%             'muscle' or 'adipose';
%
% Output:   - output image: The input image but with the segmented areas
%             colored according the atlas.
 newIm = im;

% Set the gray level of the tissue
switch tissue
    case 'bone'
        newIm(mask == 1) = 255;
    case 'prostate'
        newIm(mask == 1) = 204;
    case 'muscle'
        newIm(mask == 1) = 163;
    case 'adipose'
        newIm(mask == 1) = 63;   
    otherwise
        disp('There is no tissue with that name!');
end

end