%global
function [atlasMovedBig, Vx, Vy]= moveSmallToBig(image, atlas, atlasBig)
	 
  %create grid for small image
  sx = length(image);
  sy = length(image);
  [x, y] = meshgrid(-(sx-1)/2:(sx-1)/2, -(sy-1)/2:(sy-1)/2);
  %--------------------------------------------------------------------------
  %affine registration on the small image
  [vxim, vyim] = registerAtlas(image, atlas);
  %--------------------------------------------------------------------------
  %create grid for big image
  sx2 = length(atlasBig);
  sy2 = length(atlasBig);
  [x2, y2] = meshgrid(-(sx-(sx/sx2))/2:(sx/sx2):(sx-(sx/sx2))/2, -(sy-(sx/sx2))/2:(sy/sy2):(sy-(sx/sx2))/2);
  
  %interpolate the the small image's vector matrix to a vector matrix of the 
  %same size as the large image
  Vx = interp2(x, y, vxim, x2, y2);
  Vy = interp2(x, y, vyim, x2, y2);
  Vy(isnan(Vy)) = 0;
  Vx(isnan(Vx)) = 0;
  
  %simulate affine registration on the large image using the interpolated
  %vector matrix
  atlasMovedBig = interp2(x2, y2, atlasBig, x2+Vx, y2+Vy);
  atlasMovedBig(isnan(atlasMovedBig)) = 0;
end
