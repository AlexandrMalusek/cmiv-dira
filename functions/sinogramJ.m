function [P,r] = sinogramJ(I,thetavec,rvec,filter)
%SINOGRAMD Computes a Sinogram i.e. the Radon transform.
%   The SINOGRAMD function computes the Radon transform, which is the
%   projection of the image intensity along a radial line
%   oriented at a specific angle.
%
%   R = SINOGRAMD(I,THETA) returns the Radon transform of the
%   intensity image I for the angle THETA degrees. If THETA is a
%   scalar, the result R is a column vector containing the Radon
%   transform for THETA degrees. If THETA is a vector, then R is
%   a matrix in which each column is the Radon transform for one
%   of the angles in THETA. 
%
%   R = SINOGRAMD(I,THETA,FILTER) returns a Radon transform with the
%   detector distance equal to the pixel distance.
%   FILTER is the interpolation filter (for details, see the code):
%   1 = nearest neighbour
%   2 = linear interpolation
%   3 = cubic interpolation
%   4 = sincot, h=0.77
%   5 = sincot, h=0.6
%
%   The number of points the projection is computed as:
%   MNmax = max(size(I));
%   rLast = ceil(sqrt(2)*(MNmax-1)/2)*2+1 + 3;
%
%   This number is sufficient to compute the projection at unit
%   intervals, even along the diagonal.
%
%   [R,Xp] = SINOGRAMD(...) returns a vector Xp containing the radial
%   coordinates corresponding to each row of R.
%
%   Remarks
%   -------
%   The radial coordinates returned in Xp are the values along
%   the x-prime axis, which is oriented at THETA degrees
%   counterclockwise from the x-axis. The origin of both axes is
%   the center pixel of the image, which is defined as:
%
%        floor((size(I)+1)/2)
%
%   For example, in a 20-by-30 image, the center pixel is
%   (10,15).
%
%   Example
%   -------
%       I = zeros(100,100);
%       I(25:75, 25:75) = 1;
%       theta = 0:180;
%       [R,xp] = sinogramJD(I,theta,4);
%       imshow(theta,xp,R,[],'n')
%       xlabel('\theta (degrees)')
%       ylabel('x''')
%       colormap(hot), colorbar
%
% The routine is based on Matlabs iradon.m.                  
% Written by Maria Magnusson Seger 2003-02
% Updated by Maria Magnusson Seger 2003-03
%  

%error(nargchk(3,3,nargin))

  % 0 = Matlab, 1 = C, 2 = OpenMP
  global useCode
  switch (useCode)
    case 2
      [P,r] = openmp_sinogramJc(double(I),thetavec,rvec,filter);
      return;
  end

[P,r] = sinogramJc(double(I),thetavec,rvec,filter);
