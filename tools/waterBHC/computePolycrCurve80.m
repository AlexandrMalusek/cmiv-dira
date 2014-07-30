%=======================================================================
% Compute polychromatic and monocromatic curves for water beamhardening.
% (Based on a program by Arif Muhammad.)
% Maria 2012-02-14
% Updated by Maria 2013-10-25 (larger WBHC tables)
% Updated by Maria 2014-07-22 (drasim spectra)
%=======================================================================

% Read spectra
%-------------
load spectra_slice113B.mat
spect = zeros(76,2);
spect(:,1) = currSpectLow(1:76,1);
spect(:,2) = currSpectLow(1:76,2);

% Read water attenuation data
%----------------------------
Wdens = 1; % water density
load water.mat

% Save to matrix
%---------------
F = zeros(80,3);
F(:,1) = 1:80;
k = 1;
pos = 1;
while abs(spect(pos,1)- k) > 0.1
  k = k + 1;
end
while (k < 81)
  F(k,2) = spect(pos,2);
  k = k + 1;
  pos = pos + 1;
end 
k = 1;
pos = 1;
while (k < 81)
  while abs(1000 * f(pos,1)- k) > 0.1
    pos = pos + 1;
  end
  F(k,3) = Wdens * (f(pos,2) + f(pos,3) + f(pos,4));
  k = k + 1;
end 

% Plot spectrum and water attenuation
%------------------------------------
figure(2)
subplot(1,3,1), plot(F(:,1),100*Wdens*F(:,3)) 
title('attenuation coefficient')
xlabel('E [keV]'); ylabel('\mu(E) [1/m]');
axis([10 80 0 100])
grid
subplot(1,3,2),plot(F(:,1),F(:,2))
title('Spectrum')
xlabel('E [keV]'), ylabel('N(E)'),
xlim([10 80]), grid

% Some initializations
%---------------------
E    = F(:,1);          % spectrum energies 0-80(kev)   
N    = F(:,2);          % number of photons
mukV = F(:,3);          % attenuation coefficient, diff kV
I    = zeros(1,100);    % initialize intensity at detector
IkV  = zeros(1,79);     % initialize intensity at detector, diff KV

% Compute incoming intensity
%---------------------------
IkV = 0*IkV;
for k = 2:79
  IkV(k) = E(k) * N(E(k)) * (E(k+1) - E(k-1));
end
I0 = sum(IkV)/2;
  
% Loop over different distances (1-100cm) and 
% compute intensity at detector
%--------------------------------------------
for dist = 1:100;
  for k = 2:79;
    IkV(k) = E(k) * N(E(k)) * exp(-mukV(E(k))*dist) * (E(k+1) - E(k-1));
  end
  I(dist) = sum(IkV)/2;
end  
  
% Compute effective attenuation coefficient 
% (according to Arif o Maria)  
%------------------------------------------
mu_E = zeros(1,79);
for k = 2:79;     
  mu_E(k) = E(k) * N(E(k)) * mukV(E(k)) * (E(k+1) - E(k-1));
end
muEff = sum(mu_E)/(2*I0);
    
% Compute and plot polychromatic and monocromatic curves
%-------------------------------------------------------
distax = 1:100;
monocr = muEff * distax;
polycr = -log(I/I0);
subplot(1,3,3); plot(distax, polycr, '*r', distax, monocr, '*b');
title('Polychromatic and Monochromatic curves');
xlabel('dist [cm]');
ylabel('line integral');
legend('polychromatic','monochromatic')
grid;

save polycr80 polycr
save muEff80 muEff
