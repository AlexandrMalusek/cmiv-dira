% Nd: number of detector elements
% Np: number of projections
% Nr: reconstructed image size is Nr x Nr
% Ni: number of iterations
% Ncl: number of channels of Ul spectrum 
% Nch: number of channels of Uh spectrum
% Nt2: number of material doublets
% Nt3: number of material triplets

classdef PhantomModelData < handle
  properties
    numbiter      % number of iterations
    eEL           % low effective energy in keV
    eEH           % high effective energy in keV
    % Projections and reconstructed images
    projLow       % [Nd x Np double] projections for Ul
    projHigh      % [Nd x Np double] projections for Uh
    projLowBH     % [Nd x Np double] WBHC projections for Ul
    projHighBH    % [Nd x Np double] WBHC projections for Uh
    recLowSet     % {Ni x 1 cell}[Nr x Nr double]
    recHighSet    % {Ni x 1 cell}[Nr x Nr double]
    % Two-material decomposition (MD2) data
    p2MD          % boolean, perform MD2?
    name2         % {Nt2 x 1 cell}
    Dens2         % {Nt2 x 1 cell}[1 x 2 double] tabulated mass density
    Att2          % {Nt2 x 1 cell}[2 x 2 double] tabulated LACs
    tissueOrder2  %
    densSet       % {Ni x 1 cell}
    Wei2Set       % {Ni x 1 cell}[Nr x Nr double] calculated mass fractions
    % Three-material decomposition (MD3) data
    p3MD          % boolean, perform MD3?
    name3         % {Nt3 x 1 cell}
    Dens3         % {Nt3 x 1 cell}[1 x 3 double] tabulated mass density
    Att3          % {Nt3 x 1 cell}[2 x 3 double] tabulated LACs
    tissueOrder3  %
    Wei3Set       % {Ni x 1 cell}[Nr x Nr double] calculated mass fractions
    % Linear attenuation coefficients for MD2 and MD3
    muLow         % [Ncl x (Nt2+Nt3) double] LACs of doublets and triplets at spectrum energies
    muHigh        % [Nch x (Nt2+Nt3) double] LACs of doublets and triplets at spectrum energies
    % Post processing three-material decomposition data
    name3SA
    Dens3SA
    Att3SA
    mu3LowSA
    mu3HighSA
    maskSA
    Wei3SA
    WeiAv
  end

  methods
    function PlotMu(pmd, varargin)
      % Plot energy spectra
      semilogy(pmd.muHigh, '+-', varargin{:})
      title('LAC')
      xlabel('energy channel number')
      ylabel('LAC (1/m)')
    end

    function PlotRecLacImages(pmd, iter)
      % Plot reconstructed images (maps of LACs)
      %
      % iter: iteration number (0, ...)

      plotRecLacImages(pmd.recLowSet{iter+1}, pmd.eEL, pmd.recHighSet{iter+1}, pmd.eEH, iter);
    end

    function PlotMassFractionsFromMd3(pmd, iter)
      % Plot mass fractions from three-material decomposition
      %
      % iter: iteration number (0, ...)

      % Define a colormap
      jettmp = colormap(jet(128+32));
      jetmod = jettmp(1+16:128+16,:);
      jetmod(1,:) = jettmp(1,:); 
      jetmod(128,:) = jettmp(128+32,:); 

      for i = 1:length(pmd.Wei3Set{iter+1})
	figure()
        hold off;
        colormap(jetmod);
	subplot(1,3,1), imagesc(100 * pmd.Wei3Set{iter+1}{i}(:, :, 1), [-50 150]);
	axis image; axis off; colorbar('SouthOutside'); title(pmd.name3{i}{1});
	subplot(1,3,2), imagesc(100 * pmd.Wei3Set{iter+1}{i}(:, :, 2), [-50 150]);
	axis image; axis off; colorbar('SouthOutside'); title(pmd.name3{i}{2});
	subplot(1,3,3), imagesc(100 * pmd.Wei3Set{iter+1}{i}(:, :, 3), [-50 150]);
	axis image; axis off; colorbar('SouthOutside'); title(pmd.name3{i}{3});
      end
    end
    
    function PlotMassFractionsFromMd2(pmd, iter)
      % Plot mass fractions from two-material decomposition and the mass density
      %
      % iter: iteration number (0, ...)

      % Define a colormap
      jettmp = colormap(jet(128+32));
      jetmod = jettmp(1+16:128+16,:);
      jetmod(1,:) = jettmp(1,:); 
      jetmod(128,:) = jettmp(128+32,:); 

      for i = 1:length(pmd.Wei2Set{iter+1})
	figure()
        hold off
	colormap(jetmod);
	subplot(1,3,1), imagesc(100 * pmd.Wei2Set{iter+1}{i}(:, :, 1), [-50 150]);
	axis image; axis off; colorbar('SouthOutside'); title(pmd.name2{i}{1});
	subplot(1,3,2), imagesc(100 * pmd.Wei2Set{iter+1}{i}(:, :, 2), [-50 150]);
	axis image; axis off; colorbar('SouthOutside'); title(pmd.name2{i}{2});
	subplot(1,3,3), imagesc(pmd.densSet{iter+1}{1});
	axis image; axis off; colorbar('SouthOutside'); title('density (g/cm^3)');
      end
    end

    function [meanReg, stdReg] = meanAndStdInRegions(pmd, centerX, centerY, radius, color, iter)
      % meanAndStdInRegions

      N0 = size(pmd.recLowSet{iter+1}, 1);  % image size of pixels
      r2 = radius.^2;                       % radius squared

      % Process each region separately
      for i = 1:length(radius)
        [ix,iy] = meshgrid(1:N0, 1:N0);
        R2 = (ix - centerX(i)).^2 + (iy - centerY(i)).^2;
	% Evaluate regions reconstructed at Ul
        v = pmd.recLowSet{iter+1}(R2 < r2(i));
        meanReg(i, 1) = mean(v);
        stdReg(i, 1) = std(v);
        % Evaluate regions reconstructed at Uh
        v = pmd.recHighSet{iter+1}(R2 < r2(i));
        meanReg(i, 2) = mean(v);
        stdReg(i, 2) = std(v);
      end

      % Plot reconstructed images and regions
      figure();

      % Image for Ul
      subplot(1,2,1), imagesc(pmd.recLowSet{iter+1}, [14 30]);
      hold on;
      axis image; axis off; colorbar('horiz');
      title(sprintf('LAC (1/m), E=%.1f keV, Ni=%d', pmd.eEL, iter), 'fontsize', 11);
      colormap('gray');
      t = 0:pi/50:2*pi;
      for i = 1:length(radius)
        plot(radius(i)*sin(t)+centerX(i), radius(i)*cos(t)+centerY(i), color{i});
      end
      hold off;
      
      % Image for Uh
      subplot(1,2,2), imagesc(pmd.recHighSet{iter+1}, [14 30]);
      hold on;
      axis image; axis off; colorbar('horiz');
      title(sprintf('LAC (1/m), E=%.1f keV, Ni=%d', pmd.eEH, iter), 'fontsize', 11);
      colormap('gray');
      for i = 1:length(radius)
        plot(radius(i)*sin(t)+centerX(i), radius(i)*cos(t)+centerY(i), color{i});
      end
      hold off;
    end

    function cn = GetCondNumMdL3(pmd)
      % Return condition number of MD3 using LACs to get volume fractions.
      % Compute the condition number for each triplet.

      for i = 1:length(pmd.Att3)
	Av = zeros(3,3);           % System matrix for volume fractions
	Av(1:2,:) = pmd.Att3{i};
	Av(3,:) = [1 1 1];
	cn(i) = cond(Av);
	fprintf('Condition number for triplet (%s, %s, %s) = %g\n',...
          char(pmd.name3{i}(1)), char(pmd.name3{i}(2)), char(pmd.name3{i}(3)), cn(i));
      end
    end

    function cn = GetCondNumMdM3(pmd)
      % Return condition number of MD3 using MACs to get mass fractions.
      % Compute the condition number for each triplet.

      for i = 1:length(pmd.Att3)
	Am = zeros(3,3);           % System matrix for mass fractions
	Am(1,:) = pmd.Att3{i}(1,:) ./ pmd.Dens3{i};
	Am(2,:) = pmd.Att3{i}(2,:) ./ pmd.Dens3{i};
	Am(3,:) = [1 1 1];
	cn(i) = cond(Am);
	fprintf('Condition number for triplet (%s, %s, %s) = %g\n',...
          char(pmd.name3{i}(1)), char(pmd.name3{i}(2)), char(pmd.name3{i}(3)), cn(i));
      end
    end

  end % methods

end % classdef
