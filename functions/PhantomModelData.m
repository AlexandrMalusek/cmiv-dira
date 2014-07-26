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
    projLow       % [Nd x Np double] projections for Ul
    projHigh      % [Nd x Np double] projections for Uh
    projLowBH     % [Nd x Np double] WBHC projections for Ul
    projHighBH    % [Nd x Np double] WBHC projections for Uh
    name2         % {Nt2 x 1 cell}
    name3         % {Nt3 x 1 cell}
    Dens2         % {Nt2 x 1 cell}[1 x 2 double] mass density for MD2
    Att2          % {Nt2 x 1 cell}[2 x 2 double] LAC for MD2
    Dens3         % {Nt3 x 1 cell}[1 x 3 double] mass density for MD3
    Att3          % {Nt3 x 1 cell}[2 x 3 double] LAC for MD3
    muLow         % [Ncl x (Nt2+Nt3) double] LACs of doublets and triplets at spectrum energies
    muHigh        % [Nch x (Nt2+Nt3) double] LACs of doublets and triplets at spectrum energies
    tissueOrder2  %
    tissueOrder3  %
    recLowSet     % {Ni x 1 cell}[Nr x Nr double]
    recHighSet    % {Ni x 1 cell}[Nr x Nr double]
    densSet       % {Ni x 1 cell}
    Wei2Set       % {Ni x 1 cell}
    Wei3Set       % {Ni x 1 cell}
  end

  methods
    function PlotMu(pmd, varargin)
      % Plot energy spectra
      semilogy(pmd.muHigh, '+-', varargin{:})
      title('LAC')
      xlabel('energy channel number')
      ylabel('LAC (1/m)')
    end

    function PlotRecLac(pmd, iter)
      % Plot reconstructed images (maps of LACs)
      %
      % iter: iteration number (0, ...)

      % Scaling for plottning
      minLow = 17;
      maxLow = 33;
      minHigh = 15;
      maxHigh = 23;
  
      figure();
      subplot(1,2,1);
      imagesc(pmd.recLowSet{iter+1}, [minLow maxLow]);
      axis image; colorbar('horiz'); axis off;
      title(sprintf('80kV, recon. No: %d', iter));
      colormap('bone');
      
      subplot(1,2,2);
      imagesc(pmd.recHighSet{iter+1}, [minHigh maxHigh]);
      axis image; colorbar('horiz'); axis off;
      title(sprintf('Sn140kV, recon. No: %d', iter));
      colormap('bone');
    end

    function PlotMassFractionsFromMd3(pmd, iter)
      % Plot mass fractions from three-material decomposition
      %
      % iter: iteration number (0, ...)

      for i = 1:length(pmd.Wei3Set{iter+1})
	figure()
	subplot(1,3,1), imagesc(100 * pmd.Wei3Set{iter+1}{i}(:, :, 1), [0 100]);
	axis image; axis off; colorbar('SouthOutside'); title(pmd.name3{i}{1});
	subplot(1,3,2), imagesc(100 * pmd.Wei3Set{iter+1}{i}(:, :, 2), [0 100]);
	axis image; axis off; colorbar('SouthOutside'); title(pmd.name3{i}{2});
	subplot(1,3,3), imagesc(100 * pmd.Wei3Set{iter+1}{i}(:, :, 3), [0 100]);
	axis image; axis off; colorbar('SouthOutside'); title(pmd.name3{i}{3});
      end
    end
    
    function PlotMassFractionsFromMd2(pmd, iter)
      % Plot mass fractions from two-material decomposition and the mass density
      %
      % iter: iteration number (0, ...)

      for i = 1:length(pmd.Wei2Set{iter+1})
	figure()
	subplot(1,3,1), imagesc(100 * pmd.Wei2Set{iter+1}{i}(:, :, 1), [0 100]);
	axis image; axis off; colorbar('SouthOutside'); title(pmd.name2{i}{1});
	subplot(1,3,2), imagesc(100 * pmd.Wei2Set{iter+1}{i}(:, :, 2), [0 100]);
	axis image; axis off; colorbar('SouthOutside'); title(pmd.name2{i}{2});
	subplot(1,3,3), imagesc(pmd.densSet{iter+1}{1});
	axis image; axis off; colorbar('SouthOutside'); title('density (g/cm^3)');
      end
    end

  end % methods

end % classdef
