function q = ramp(g, dalfa)
  % q = ramp(g, dalfa)
 
  % rampfilter the columns of g
  %  g:       input to be filtered
  %  dalfa:   angle between the rays
  %
  %  w = cos(pi*(-M2:M2)'/(2*M2)).^2; % low pass weighting

  [M,N] = size(g);

  M2 = 2^ceil(log2(M));
  r = zeros(2*M2,1);
  if nargin==1
    r(M2+1) = 1/4;
    r(2:2:end) = -1 ./(pi^2*(-M2+1:2:M2-1).^2);
  end
  if nargin==2   
    r(M2+1) = 1/(8*dalfa^2);
    r(2:2:end) = -1./(2*pi^2*sin((-M2+1:2:M2-1)*dalfa).^2);
  end
  
  w = cos(pi*(-M2:M2)'/(2*M2)).* cos(pi*(-M2:M2)'/(2*M2)); % low pass weighting
    R = (w(1:end-1) .* fftshift(fft(ifftshift(r))))*ones(1,N);
    
  q = zpadcol(g,2*M2);
  q = real(fftshift(ifft(ifftshift (R.*fftshift(fft(ifftshift(q)))))));
  q = zpadcol(q,M);
end

function y = zpadcol(x,m2)
  %
  % y = zpadcol(x,m2)
  % pad or unpad
  %
  
  [m,n] = size(x);

  if mod(m-m2,2) == 0
    offset = abs((m-m2)/2);
  else
    offset = (abs(m-m2)+1)/2;
  end
  
  if m2>m
    y = [zeros(offset,n); x; zeros(m2-m-offset,n)];
  else
    y = x(offset+1:offset+m2,:);
  end
end
