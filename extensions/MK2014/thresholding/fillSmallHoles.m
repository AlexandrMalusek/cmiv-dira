%fill small holes
function image = fillSmallHoles(original, size)
	 
  filled = imfill(original, 'holes');
  holes = filled & ~original;
  bigholes = bwareaopen(holes, size);
  smallholes = holes & ~bigholes;
  new = original | smallholes;
  image = new;
end
