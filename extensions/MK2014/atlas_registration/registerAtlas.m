function [vxim, vyim] = registerAtlas(image, atlas)
	 
  %load quadrature filters(2)
  load('quadratureFilters.mat')
  
  %create a grid for the image
  [sy, sx] = size(image);
  [x, y] = meshgrid(-(sx-1)/2:(sx-1)/2, -(sy-1)/2:(sy-1)/2);
  
  ptot = zeros(6, 1);
  
  %create the vector matrixes
  [sy2, sx2] = size(image);
  vxim = zeros(sy2, sx2);
  vyim = zeros(sy2, sx2);
  
  %loop until difference is small enough
  diff = 100;
  loops = 0;
  while diff >= 0.002 || loops == 200
	
    %move the atlas to closer resemble the input image
    atlasMoved = interp2(x, y, atlas, x+vxim, y+vyim);
    atlasMoved(isnan(atlasMoved)) = 0;
    
    %create constants from the quadrature filters
    for k = 1:2
      q1{k} = conv2(image, f{k}.v, 'same');
      q2{k} = conv2(atlasMoved, f{k}.v, 'same');
      
      phi1{k} = angle(q1{k});
      phi2{k} = angle(q2{k});
      
      deltaT{k} = angle(q2{k} .* conj(q1{k}));
      c{k} = sqrt(abs(q1{k} .* q2{k})) .* cos(deltaT{k}/2)^2;
    end
    
    %create matrix used for the calculation of the registration
    A = [sum(c{1}(:).^2) 0 sum(x(:).*c{1}(:).^2) sum(y(:).*c{1}(:).^2) 0 0;...
        0 sum(c{2}(:).^2) 0 0 sum(x(:).*c{2}(:).^2) sum(y(:).*c{2}(:).^2);...
        0 0 sum(x(:).^2.*c{1}(:).^2) sum(x(:).*y(:).*c{1}(:).^2) 0 0;...
        0 0 0 sum(y(:).^2.*c{1}(:).^2) 0 0;...
        0 0 0 0 sum(x(:).^2.*c{2}(:).^2) sum(x(:).*y(:).*c{2}(:).^2);...
        0 0 0 0 0 sum(y(:).^2.*c{2}(:).^2)];
    
    Dx = phi1{1} - phi2{1};
    Dy = phi1{2} - phi2{2};
    
    Atot = A + triu(A, 1)';
    
    h = [sum(c{1}(:).^2.*Dx(:)); sum(c{2}(:).^2.*Dy(:)); sum(x(:).*c{1}(:).^2.*Dx(:));...
      sum(y(:).*c{1}(:).^2.*Dx(:)); sum(x(:).*c{2}(:).^2.*Dy(:)); sum(y(:).*c{2}(:).^2.*Dy(:))];
    
    p = Atot \ h;
    ptot = ptot + p;
    vx = 1*ptot(1) + x*ptot(3) + y*ptot(4);
    vy = 1*ptot(2) + x*ptot(5) + y*ptot(6);
    
    vxim = reshape(vx, sx, sy);
    vyim = reshape(vy, sx, sy);
    
    diff = sum(abs(p));
    loops = loops + 1;
  end
end
