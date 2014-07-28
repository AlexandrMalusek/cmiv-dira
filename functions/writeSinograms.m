%% Write sinograms

% Read drasim generated projections and perform:
% - parallel rebinning
% - water beam hardening correction
% - image reconstruction (for display only)
% Save resulting sinograms and drasim generated spectra
%
% (Based on a program by Arif Muhammad.)
% Maria 2012-02-22
% Oscar Grandell 2012-08-27 - Heavy changes
% Robin Westin 2012-10-18 - Updated
% Maria Magnusson 2013-07-17 - Updated
% Alexandr Malusek 2014-07-13 - Simplified

disp('Rebinning ...');

for i = 1:2
  % Process data for low and high tube voltages

  if i == 1
    spect = 'Low';
  else
    spect = 'High';
  end
  % load polycromatic curve and effective attenuation coefficient
  %--------------------------------------------------------------
  if strcmp(spect, 'Low');
    load polycr80
    load muEff80
    minv = 15;  % 1/m
    maxv = 55;  % 1/m
  end
  if strcmp(spect, 'High');
    load polycr140Sn
    load muEff140Sn
    minv = 15;  % 1/m
    maxv = 35;  % 1/m
  end

  % Using drasim data
  %------------------
  if strcmp(spect, 'Low');
    drasim = flipud(create_sinogram(0, smd.M0, smd.N0)'); % Create sinogram of data
  elseif strcmp(spect, 'High');
    drasim = flipud(create_sinogram(1, smd.M0, smd.N0)'); % Create sinogram of data
  end
  drasim = -log(drasim/max(max(drasim)));
  
  % Rebinning of fan beam data
  %---------------------------
  rebsim = rebinning(drasim, smd.L, smd.dt0, smd.dfi0, 1, 0, smd.dt1, smd.dfi1, smd.N1, smd.M1);
  
  % Projection data (save for iterative loop)
  if strcmp(spect, 'Low');
    projLow = rebsim; 
  elseif strcmp(spect, 'High');
    projHigh  = rebsim;
  end
  recPolysim = parFB(rebsim, smd.dt1, smd.gamma); % Parallel filtered BP
  
  %  Plots of data without BHC
  %---------------------------
  figure()
  subplot(221); imagesc(drasim); title('Fan beam sinogram'); colorbar;
  subplot(222); imagesc(rebsim); title('Parallel beam sinogram'); colorbar;
  subplot(223); imagesc(recPolysim); title('Reconstructed LAC (1/m)'); colorbar;
  subplot(224); plot(recPolysim(smd.N0/2,:)); title('LAC profile (1/m)');
  axis([1 smd.N1 minv maxv]);
  grid on;
  
  % Using CT-numbers
  %-----------------
  recPolyHU = 1000 * (recPolysim  - muEff*100)/(muEff*100);
  figure(); imagesc(recPolyHU);title('Reconstructed CT numbers (HU)'); colorbar
  
  % Preforming BHC
  %---------------
  dist = zeros(smd.N1, smd.M1);      
  for  n = 1:smd.M1
    for  m = 1:smd.N1
      value = rebsim(m,n);
      for k = 1:99
        if (value < polycr(k) && k==1)
          dist(m,n) = value / polycr(1); 
          break;    
        else
          if (value>=polycr(k)) && (value<polycr(k+1))
            dist(m,n) = k + (value-polycr(k)) / (polycr(k+1)-polycr(k)); 
            break;
          end
        end
      end
    end
  end 
  
  BHcorr = muEff * dist ./ (rebsim+eps); 
  rebprojBHcorr = rebsim .* BHcorr;
  
  % Beam hardening corrected projection data (save for iterative loop)
  %-------------------------------------------------------------------
  if strcmp(spect, 'Low');
    projLowBH = rebprojBHcorr; 
  elseif strcmp(spect, 'High');
    projHighBH  = rebprojBHcorr;
  end
  
  fId = fopen(sprintf('spektren/spk1_0.%d.spc', i-1), 'r');
  spectCellTmp = textscan(fId, '%f %f', [2 inf]);
  if strcmp(spect, 'Low');
    currSpectLow(:, 1) = spectCellTmp{1};
    currSpectLow(:, 2) = spectCellTmp{2};
  elseif strcmp(spect, 'High');
    currSpectHigh(:, 1) = spectCellTmp{1};
    currSpectHigh(:, 2) = spectCellTmp{2};
  end
  fclose(fId);
end

save(sinogramsFileName, 'projLow', 'projHigh');
save(sinogramsBhFileName, 'projLowBH', 'projHighBH');
save(spectraFileName, 'currSpectLow', 'currSpectHigh');
