function q = rampwindowMtf(p, rVec, n)
  % rampwindowMtf
  %
  %===============================================================
  % Window function MTF to be used together with ramp-filter.
  % The Window function is divided by cos^n, (=0 gives no effect)
  % Maria Magnusson, 2014-07-08
  % Alexander Örtenberg 2014-11-12
  %================================================================

  % variables
  %----------
  Nr = length(rVec);
  delta = rVec(2)-rVec(1);
  Raxis = (-Nr:Nr-1)'/(2*Nr*delta);
  Nphi = size(p,2);
  odd = rem(Nr,2);
  
  % design a MTF weighting window
  %----------------- -------------
  MTF = load('MTF.dat', '-ascii');
  W = 0*Raxis;
  for k = 1:2*Nr
    W(k) = mtfWindow(abs(Raxis(k)), MTF);
  end;
  
  % divide by cos^n, (=0 gives no effect)
  %--------------------------------------
  W = W./cos(pi*Raxis*delta/2).^n;
  
  % figure(77)
  % plot(Raxis,W,'r')
  % axis([-1/(2*delta),1/(2*delta),0,1.1])
  % grid
  
  % Apply the window in the Fourier domain
  %---------------------------------------
  Wmat = W  * ones(1,Nphi);
  pp = [zeros((Nr+odd)/2, Nphi); p; zeros((Nr-odd)/2, Nphi)];
  P = fftshift(fft(ifftshift(pp), [], 1));
  Q = P .* Wmat;
  q = real(fftshift(ifft(ifftshift(Q), [], 1)));
  q = q((Nr+odd)/2+1:(Nr+odd)/2+Nr, :);
end
