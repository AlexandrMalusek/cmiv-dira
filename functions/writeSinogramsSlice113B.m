%===================================================================
% Projection data and spectrum are loaded from drasim with 
% real geometry including quarter offset.
% Note that drasim data is swapped here.
% Rebinning and reconstruction is performed.
% Water beam hardening correction (WBHC) followed by reconstruction
% is also performed.
% (Based on a program by Arif Muhammad.)
% Maria Magnusson 2012-02-22
% Oscar Grandell 2012-08-27 - Heavy changes
% Robin Westin 2012-10-18 - Updated
% Maria Magnusson 2014-06-11 - Changes to real geometry, better WBHC
% Maria Magnusson 2014-07-14 - Changes to cubic interpolation during 
%                              rebinning, better plots
%====================================================================

% Initialization of projection vectors
% ====================================
detfanVec = (-(smd.N0-1)/2+1.25:(smd.N0-1)/2+1.25)*smd.dt0; % for swapped detector
betaVec    = -1 * (0:smd.M0-1) * smd.dfi0; % start from 0 going negative...
detparVec  = (-(smd.N1-1)/2:(smd.N1-1)/2)*smd.dt1;
phiVec     = -1 * (0:smd.M1-1) * 180/smd.M1 - smd.gamma*180/pi;

% 90 deg extra for simpler rebinning 
% ----------------------------------
betaVecEXTRA = -1 * (0:(smd.M0*1.25)-1) * smd.dfi0; 

% 180 deg extra for simpler rebinning 
% -----------------------------------
phiVecEXTRA  = -1 * (0:2*smd.M1-1) * 180/smd.M1 - smd.gamma*180/pi; 

% Initialization of mask
% ======================
[x,y] = meshgrid(-(smd.N1-1)/2:(smd.N1-1)/2,-(smd.N1-1)/2:(smd.N1-1)/2);
mask = (x.^2 + y.^2) <= ((smd.N1-1)/2)^2;

% Perform the procedure for both 'Low' and 'High' spectral data
% =============================================================
for i = 1:2
    
  if i == 1
    spect = 'Low';
  else
    spect = 'High';
  end
  
  % Load spectral data and effective attenuation coefficient
  % --------------------------------------------------------
  if strcmp(spect, 'Low');
    load polycr80
    load muEff80
    minv = 19;
    maxv = 25;
    maxvHU = 3000;
  end
  if strcmp(spect, 'High');
    load polycr140Sn
    load muEff140Sn
    minv = 16;
    maxv = 22;
    maxvHU = 2000;
  end
  
  % Load drasim data
  % ----------------
  if strcmp(spect, 'Low');
    drasimOrig = create_sinogram(0,smd.M0,smd.N0)'; % Create sinogram of data
  elseif strcmp(spect, 'High');
    drasimOrig = create_sinogram(1,smd.M0,smd.N0)'; % Create sinogram of data
  end
  drasimOrig = -log(drasimOrig/max(max(drasimOrig)));
  
  if i==1
    figure(1)
    imagesc(drasimOrig);
    title('Original drasim Low energy fanbeam sinogram');
    colorbar;
    xlabel('projection angle no');
    ylabel('detector element no');
  end
  
  % Swap drasim data in the detector direction
  % ------------------------------------------
  drasimOrig = flipud(drasimOrig);
  
  % Extend fanbeam sinogram with 90 deg
  % -----------------------------------
  drasim = [drasimOrig, drasimOrig(:,1:smd.M0/4)];
  
  % Rebinning of fanbeam data
  % -------------------------
  [F,t]   = meshgrid(betaVecEXTRA, detfanVec);
  [Fn,tn] = meshgrid(phiVecEXTRA,  detparVec);
  rebsim = interp2(F, t, drasim, Fn-180*asin(tn/smd.L)/pi, asin(tn/smd.L), 'cubic');
  loc = find(isnan(rebsim)==1);   % find NaN:s and replace them with 0
  rebsim(loc) = 0;
  
  % Merge sinogram from 360 to 180 deg
  % ----------------------------------
  rebsim = 0.5 * (rebsim(:,1:smd.M1) + rebsim(end:-1:1,smd.M1+1:2*smd.M1));
  
  % Copy projection data for later use in DIRA
  % ------------------------------------------
  if strcmp(spect, 'Low');
    projLow = rebsim; 
  elseif strcmp(spect, 'High');
    projHigh  = rebsim;
  end
  
  % Perform parallel FBP with mask
  % ------------------------------
  if 1==2
    rebsimfilt = rampwindowb(rebsim, detparVec, 'Cosine', 4, 1.8);
  else
    n=4;
    rebsimfilt = rampwindowMtf(rebsim, detparVec*100, n);
  end
  recPolysim = iradon(rebsimfilt, phiVec, 'linear', 'Ram-Lak', 1, smd.N1)/smd.dt1;  
  recPolysim = recPolysim .* mask;

  % Plot data without WBHC
  %-----------------------
  if i==1
    figure(12)
    subplot(221); imagesc(drasim);
    title('drasim data, swapped and extended');colorbar;
    subplot(222); imagesc(rebsim);
    title('rebinned sinogram');colorbar;
    subplot(223); imagesc(recPolysim); axis image
    title('Low energy reconstruction [LAC]');colorbar;
    subplot(224); plot(recPolysim((smd.N1+1)/2,:));
    title('reconstruction, central cut [LAC]');
    axis([1 smd.N1 minv maxv]); grid on;
  elseif i==2
    figure(22)
    subplot(221); imagesc(drasim);
    title('drasim data, swapped and extended');colorbar;
    subplot(222); imagesc(rebsim);
    title('rebinned sinogram');colorbar;
    subplot(223); imagesc(recPolysim); axis image
    title('High energy reconstruction [LAC]');colorbar;
    subplot(224); plot(recPolysim((smd.N1+1)/2,:));
    title('reconstruction, central cut [LAC]');
    axis([1 smd.N1 minv maxv]); grid on;
  end
  
  % Perform WBHC (obs! now up to k=99)
  %-----------------------------------
  dist = zeros(smd.N1,smd.M1);      
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
  
  BHcorr = muEff * dist ./ (rebsim + eps); 
  rebprojBHcorr = rebsim .* BHcorr;
  
  % Perform parallel FBP with mask
  % ------------------------------
  if 1==2
    rebprojBHcorrfilt = rampwindowb(rebprojBHcorr, detparVec, 'Cosine', 4, 1.8);
  else
    n=4;
    rebprojBHcorrfilt = rampwindowMtf(rebprojBHcorr, detparVec*100, n);
  end
  recPolysim = iradon(rebprojBHcorrfilt, phiVec, 'linear', 'Ram-Lak', 1, smd.N1)/smd.dt1;
  recPolysim = recPolysim .* mask;
  
  % Plot data with WBHC
  % -------------------
  if i==1
    figure(13)
    subplot(222); imagesc(rebprojBHcorr);
    title('rebinned sinogram + WBHC');colorbar;
    subplot(223); imagesc(recPolysim); axis image
    title('WBHC Low energy reconstruction [LAC]');colorbar;
    subplot(224); plot(recPolysim((smd.N1+1)/2,:));
    title('reconstruction, central cut [LAC]');
    axis([1 smd.N1 minv maxv]); grid on;
    figure(14); colormap gray
    hold off
    imagesc(recPolysim);title('WBHC Low energy reconstruction [LAC]');
    axis image
    colorbar;
    caxis([minv maxv])
  elseif i==2
    figure(23)
    subplot(222); imagesc(rebprojBHcorr);
    title('rebinned sinogram + WBHC');colorbar;
    subplot(223); imagesc(recPolysim); axis image
    title('WBHC High energy reconstruction [LAC]');colorbar;
    subplot(224); plot(recPolysim((smd.N1+1)/2,:));
    title('reconstruction, central cut [LAC]');
    axis([1 smd.N1 minv maxv]); grid on;
    figure(24); colormap gray
    hold off
    imagesc(recPolysim);title('WBHC High energy reconstruction [HU]');
    axis image
    colorbar;
    caxis([minv maxv])
  end
  
  % Use HU-numbers and measure mean and std
  % ---------------------------------------
  recPolyHU = 1000 * (recPolysim  - muEff*100)/(muEff*100);
  if i==1
    figure(15); colormap gray
    hold off
    imagesc(recPolyHU);title('WBHC Low energy reconstruction [HU]');
    axis image
    colorbar;
    caxis([-100 100])
  elseif i==2
    figure(25); colormap gray
    hold off
    imagesc(recPolyHU);title('WBHC High energy reconstruction [HU]');
    axis image
    colorbar;
    caxis([-100 100])
  end   
  
  % Copy WBHC projection data for later use in DIRA
  % -----------------------------------------------
  if strcmp(spect, 'Low');
    projLowBH = rebprojBHcorr; 
  elseif strcmp(spect, 'High');
    projHighBH  = rebprojBHcorr;
  end
  
  % Copy Spectra for later use in DIRA
  % ----------------------------------
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
