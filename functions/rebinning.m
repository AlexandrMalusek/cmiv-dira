function y = rebinning(x, L, dt, dfi, fb, rot, dtn, dfin, Mn, Nn)
% function y = rebin(x, L, dt, dfi, rot, dtn, dfin, Mn, Nn)
% L: distance to source
% dt: detector element size
% dfi: angle increment degrees
% fb: clockwise or counterclockwise rotation?
%     try either fb = 0 and fb = 1
% rot: rotation center pixel nr
% dtn: new detector element size
% dfin: new angle increment
% Mn: new nr of pixels per proj
% Nn: new nr of projs

[M,N] = size(x);    % M = nr of pixels per projection
                    % N = nr of projections

if nargin < 10
  Nn = N;
end
if nargin < 9
  Mn = M;
end
if nargin < 8
  dtn = dt;
end
if nargin < 7
  dfin = dfi;
end
if nargin < 6
  rot = 0;
end

gamma = 180*atan(M/2*dt/L)/pi;		% fanbeam angle

% in- and out-grid
%-----------------
if fb == 0
  [F,t]   = meshgrid((0:N-1)*dfi, (-(M-1)/2:(M-1)/2)*dt);
  [Fn,tn] = meshgrid((0:Nn-1)*dfin+gamma, (-(Mn-1)/2:(Mn-1)/2)*dtn);
else
  [F,t]   = meshgrid((0:N-1)*dfi, ((M-1)/2:-1:-(M-1)/2)*dt);
  [Fn,tn] = meshgrid((0:Nn-1)*dfin+gamma, ((Mn-1)/2:-1:-(Mn-1)/2)*dtn);
end

% rebin
%------
y=interp2(F, t+rot, x, Fn-180*asin(tn/L)/pi, asin(tn/L));

% find NaN:s and replace them with 0
%-----------------------------------
loc = find(isnan(y)==1);
y(loc) = 0;


