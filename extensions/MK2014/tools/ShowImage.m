function ShowImage(im,fig_numb)
% Function for plotting a grayscale image. Input is the image and the
% figure number. The function set the axis to image and colormap to gray
%
% Written by Julius Jeuthe, 2015 


if (1 < size(im,3))
    error('Input should be a gray scale image')
else
    figure(fig_numb)
    imagesc(im)
    axis image;
    colormap gray
    drawnow
end
end