function [atlasMoved, atlasMovedAffine] = registration(image, atlas, bones)
	 
  %--------------------------------------------------------------------------
  %AFFINE REGISTRATION-------------------------------------------------------
  %create binary image and atlas
  [imageBinary] = createBinary(image);
  [atlasBinary] = createBinary(atlas);
  
  %rescale atlas and image
  scale = 127 / length(image);
  atlasBinarySmall = imresize(atlasBinary, scale);
  imageBinarySmall = imresize(imageBinary, scale);
  
  %affine registration on the small image, translation to large image
  [~, vxim, vyim] = moveSmallToBig(imageBinarySmall, atlasBinarySmall, atlasBinary);
  
  %create a grid for the small image
  sx = length(imageBinarySmall);
  sy = length(imageBinarySmall);
  [x, y] = meshgrid(-(sx-1)/2:(sx-1)/2, -(sy-1)/2:(sy-1)/2);
  
  %create a grid for the large atlas
  sx2 = length(atlas);
  sy2 = length(atlas);
  [x2, y2] = meshgrid(-(sx-(sx/sx2))/2:(sx/sx2):(sx-(sx/sx2))/2, -(sy-(sx/sx2))/2:(sy/sy2):(sy-(sx/sx2))/2);
  
  %the affine registration translated to the atlas
  atlasMovedAffine = interp2(x2, y2, double(atlas), x2+vxim, y2+vyim);
  
  %--------------------------------------------------------------------------
  %NON-RIGID BONE REGISTRATION-----------------------------------------------
  
  [~, bonesAtlas] = histc(atlasMovedAffine, [220 256]);
    
  bonesAtlas = double(bonesAtlas);
  bones = double(bones);
  
  %create grid for the large image
  sx = length(image);
  sy = length(image);
  [x3, y3] = meshgrid(-(sx-1)/2:(sx-1)/2, -(sy-1)/2:(sy-1)/2);
  
  [vxim2, vyim2]=registerAtlas(bones, bonesAtlas);
  atlasMovedNonRigid = interp2(x3, y3, atlasMovedAffine, x3+vxim2, y3+vyim2);
  bonesAtlasMovedNonRigid = interp2(x3, y3, bonesAtlas, x3+vxim2, y3+vyim2);
  
  diceBonesAffine = dice(uint8(bones), uint8(bonesAtlas));
  diceBonesNonRigid = dice(uint8(bones), uint8(bonesAtlasMovedNonRigid));
  
  if diceBonesAffine >= diceBonesNonRigid
    atlasMoved = atlasMovedAffine;
  else
    atlasMoved = atlasMovedNonRigid;
  end
end
