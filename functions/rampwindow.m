function q = rampwindow(p, rVec, weight, n)
%=====================================================
% Window function to be used together with ramp-filter
% Maria Magnusson, 2014-07-03
%=====================================================

% variables
%----------
Nr = length(rVec);
rDelta = rVec(2) - rVec(1);
rVecDouble = (-Nr:Nr-1) * rDelta;
Nphi = size(p,2);
odd = rem(Nr,2);

% design a weighting window
%--------------------------
if weight == 'Cosine'
  W = cos(pi*(-Nr:Nr-1)'/(2*Nr));
  W = W.^n;
elseif weight == 'Nofilt'
  W = 1-0*(-Nr:Nr-1)';
else
  W = 1-0*(-Nr:Nr-1)';
end

% Apply the window in the Fourier domain
%---------------------------------------
Wmat = W  * ones(1,Nphi);
pp = [zeros((Nr+odd)/2, Nphi); p; zeros((Nr-odd)/2, Nphi)];
P = fftshift(fft(ifftshift(pp), [], 1));
Q = P .* Wmat;
q = real(fftshift(ifft(ifftshift(Q), [], 1)));
q = q((Nr+odd)/2+1:(Nr+odd)/2+Nr, :);

