%transfer curve
close all;
clear all;
im1 = imread('prostate_ct.JPG');
im1 = rgb2gray(im1);

im_d = double(dicomread('DE_Pelvis1_222'));

im2 = uint8(im_d/max(max(im_d))*255);

matched_image = histogram_matching(im_d);

% figure,imshow(im1)
% figure,imshow(im_d,[-2000,2000])
% figure,imshow(im2)
% figure,imshow(matched_image)

figure, imhist(im2)

bins = 0:255;
H = hist(im1(:), bins);
Hmod = H + eps(sum(H));
cdf_im1 = [0, cumsum(Hmod)/sum(Hmod)];

bins = 0:255;
H = hist(matched_image(:), bins);
Hmod = H + eps(sum(H));
cdf_imm = [0, cumsum(Hmod)/sum(Hmod)];

bins = 0:255;
H = hist(im2(:), bins);
Hmod = H + eps(sum(H));
cdf_im2 = [0, cumsum(Hmod)/sum(Hmod)];

figure, plot(cdf_im1);

figure, plot(cdf_imm);
hold on
plot(cdf_im2,'r');
xlabel('Pixel Value') % x-axis label
ylabel('Cumulative Probability') % y-axis label
legend('Reference Image','Input Image')

figure, plot(cdf_im2, cdf_imm);
xlabel('Input Image Pixel Value') % x-axis label
ylabel('Output Image Pixel Value') % y-axis label
