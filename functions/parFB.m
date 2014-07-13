function f = parFB(g, dtau, firstang)
%
% f = parFB(g, tau) 
% parallel filtered backprojection method
%
% g:        projection indata
% dtau:     detector element size
% firstang: angle for first projection [rad]
%
  
if (nargin == 0)|(nargin>3) 
  fprintf('1, 2 or 3 arguments is desired!\n');
  return;
end
if nargin == 1
  dtau = 1;
  firstang = 0;
end
if nargin == 2
  firstang = 0;
end

[M,N] = size(g);   % M = nr of pixels per projection, 
                   % N = nr of projections			
  
Sx = M;            % image x size
Sy = M;            % image y size
m0 = M;            % image x size in pixels
n0 = M;            % image y size in pixels
if mod(M,2)==0
  startx=-M/2;
  starty=-M/2;
else
  startx=-(M-1)/2;
  starty=-(M-1)/2;
end

% projection angles
dtheta = pi/N;
%theta = -pi/2-firstang:-dtheta:-3*pi/2+dtheta-firstang;
theta = 0+firstang:dtheta:pi-dtheta+firstang;

% detector array
if mod(M,2)==0
  det = (-M/2:M/2-1)';
else
  det = (-(M-1)/2:(M-1)/2)'; 
end

% rampfilter the projections 
q = ramp(g);

% the grid for the reconstructed image
dx = Sx/m0;
dy = Sy/n0;
[X,Y] = meshgrid(startx:dx:startx+(m0-1)*dx, starty:dy:starty+(n0-1)*dy);
f = zeros(size(X));

% backprojection
for i=1:N
  c = cos(theta(i));
  s = sin(theta(i));
  f = f + interp2(ones(M,1)*(1:N), det*ones(1,N), q, i*ones(m0,m0), X*c + Y*s);
end

% find NaN:s and replace them with 0
loc = find(isnan(f)==1);
f(loc) = 0;

f = f*dtheta/dtau;
