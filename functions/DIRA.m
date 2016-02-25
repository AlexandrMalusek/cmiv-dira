%% The iterative DIRA loop.
%
% #################################################################
% Programmed 2010 by Arif Muhammad with guidance by Maria Magnusson
% Updated 2011-04 and 2012-04 by Maria Magnusson
% Updated and 2MD+dens added 2012-08 by Oscar Grandell
% Updated 2012-10 by Robin Westin
% Updated 2013-07-17 by Maria Magnusson
% Updated 2014-07-10 by Alexandr Malusek
% #################################################################

% Get the number of iterations
pmd.savedIter = sort(pmd.savedIter);             % sort the vector
numbIter = pmd.savedIter(length(pmd.savedIter)); % get the last element 

% Check whether gDiraPlotFigures exists. If not, set it to a default va 
if (0 == exist('gDiraPlotFigures'))
  global gDiraPlotFigures;
  gDiraPlotFigures = 1;
end

%% CT scan geometry and Joseph metod
% ------------------------------------

% Rebinning of fan beam data
%---------------------------
firstang = smd.gamma;            % first angle after rebinning [rad]

% Init data for Joseph projection generation
% ------------------------------------------
Nr = smd.N1;                     % number of detector elements
Nphi = size(pmd.projLow, 2);     % number of projection angles after rebinning
degLen = 180; 		         % angular interval in degrees
pixsiz = smd.dt1;                % pixel size
Nr2 = 2*floor((Nr-1)/2)+1;       % resolution of detector elements
degVec = ((0:-1:(-Nphi+1))*(degLen/Nphi))- smd.gamma * 180/pi;
r2Vec  = (-(Nr2-1)/2:1:(Nr2-1)/2);

% Initialization of circular reconstruction area mask
% ---------------------------------------------------
[x,y] = meshgrid(r2Vec,r2Vec);
smd.mask = (x.^2 + y.^2) <= ((Nr2-1)/2)^2;

% Compute I0 for both spectra
% --------------------------
sp = zeros(1,78);
for k=2:length(smd.ELow)-1;                           
  sp(k) = (smd.ELow(k)*smd.NLow(k)*(smd.ELow(k+1)-smd.ELow(k-1)));
end
uLow = sum(sp)/2;
sp = zeros(1,138);
for k=2:length(smd.EHigh)-1;                           
  sp(k) = (smd.EHigh(k)*smd.NHigh(k)*(smd.EHigh(k+1)-smd.EHigh(k-1)));
end
uHigh = sum(sp)/2;

% Filter original projections 
pmd.projLow = rampWindowForMeasuredProjections(pmd.projLow, r2Vec);
pmd.projHigh = rampWindowForMeasuredProjections(pmd.projHigh, r2Vec);

%% Reconstruction No.0
%
fprintf('\nStarting initial reconstruction...\n')

pmd.curIterIndex = 1;
phm1 = reconstructMeasuredProjections(pmd.projLowBH, r2Vec, degVec, smd.N1, smd.dt1);
phm2 = reconstructMeasuredProjections(pmd.projHighBH, r2Vec, degVec, smd.N1, smd.dt1);

nSavedIter = length(pmd.savedIter);  % Number of saved iterations
pmd.recLowSet = cell(nSavedIter, 1);
pmd.recHighSet = cell(nSavedIter, 1);

pmd.recLowSet{pmd.curIterIndex} = phm1;
pmd.recHighSet{pmd.curIterIndex} = phm2;

% Plot reconstructed maps of linear attenuation coefficients
if gDiraPlotFigures == 1
  pmd.PlotRecLacImages(0);
  drawnow();
end

%% Inital Tissue segmentation
% 
disp('Classifying tissues...')
[tissue2 tissue3] = tissueClassification(0, smd, pmd);
pmd.tissue2Set = cell(nSavedIter);
pmd.tissue3Set = cell(nSavedIter);
pmd.tissue2Set{1} = tissue2;
pmd.tissue3Set{1} = tissue3;

% Tissue decomposition
disp('Decomposing tissues...')
AttE1mat = 0.01*phm1;		   % Change 1/m to 1/cm
AttE2mat = 0.01*phm2;

if pmd.p2MD
  dens = cell(length(tissue2), 1);
  Wei2 = cell(length(tissue2), 1);
  for i = 1:length(tissue2)
    [Wei2{i}, dens{i}] = MD2(AttE1mat, AttE2mat, pmd.Att2{pmd.tissueOrder2(i)},...
      pmd.Dens2{pmd.tissueOrder2(i)}, tissue2{i});
  end
end

if pmd.p3MD
  Wei3 = cell(length(tissue3), 1);
  for i = 1:length(tissue3)
    Wei3{i} = MD3(AttE1mat, AttE2mat, pmd.Att3{pmd.tissueOrder3(i)},...
      pmd.Dens3{pmd.tissueOrder3(i)}, tissue3{i}, 1);
  end
end

pmd.densSet = cell(nSavedIter, 1);
pmd.Wei2Set = cell(nSavedIter, 1);
pmd.Wei3Set = cell(nSavedIter, 1);
if pmd.p2MD
  pmd.densSet{pmd.curIterIndex} = dens;
  pmd.Wei2Set{pmd.curIterIndex} = Wei2;
end
if pmd.p3MD
  pmd.Wei3Set{pmd.curIterIndex} = Wei3;
end

% Plot computed mass fractions from MD2 and MD3
if gDiraPlotFigures == 1
  if pmd.p2MD
    pmd.PlotMassFractionsFromMd2(0);
  end
  if pmd.p3MD
    pmd.PlotMassFractionsFromMd3(0);
  end
  drawnow();
end

%% Iterate
%
iterno = numbIter;
for iter = 1:numbIter
  % Projection generation with Joseph
  %----------------------------------
  % For every 2 or 3 MD added this part have to be extended with a loop
  % for projection calculation, monoenergetic loop and summation and
  % polychromatic concatination.
  
  fprintf('\nStarting iteration %d...\n', iter);

  % If the previous iteration is to be saved then increment the iteration index.
  % Otherwise the current iteration will overwrite the data of the previous iteration.
  if pmd.GetIterIndex(iter-1) > 0
     pmd.curIterIndex = pmd.curIterIndex + 1;
  end

  % Calculate volume fractions v_i (Vol3) from mass fractions w_i (Wei3):
  %   v_i(x,y) = w_i(x,y) * rho(x,y) / rho_i
  %   where rho(x,y) = 1/(w_1(x,y)/rho_1 + w_2(x,y)/rho_2 + w_3(x,y)/rho_3)
  for i = 1:length(tissue3)
    if i == 1
      % Temporary workaround: Handle air outside the imaged object.
      % lipidAirMask: air outside the body AND inside the CT-scanner specific "circle"
      % tissue3{1}: the union of the body soft tissue AND the air outside the body
      % In the future a different solution will be used:
      % - air outside the body will be decomposed to air + soft tissue using 2MD
      % - CT table has to be segmented out
      lipidAirMask = ((AttE1mat < pmd.Att3{1}(1,1)) | (AttE2mat < pmd.Att3{1}(2,1))) .* tissue3{1};
      Vol3{i}(:,:,1) = Wei3{i}(:,:,1) ./ ...
        (Wei3{i}(:,:,1)/pmd.Dens3{i}(1) + Wei3{i}(:,:,2)/pmd.Dens3{i}(2) + Wei3{i}(:,:,3)/pmd.Dens3{i}(3) + eps) /...
        pmd.Dens3{1}(1) .* (tissue3{1} - lipidAirMask) + Wei3{1}(:,:,1) .* lipidAirMask;
    else
      Vol3{i}(:,:,1) = Wei3{i}(:,:,1) ./ ...
        (Wei3{i}(:,:,1)/pmd.Dens3{i}(1) + Wei3{i}(:,:,2)/pmd.Dens3{i}(2) + Wei3{i}(:,:,3)/pmd.Dens3{i}(3) + eps) /...
        pmd.Dens3{1}(1); 
    end
    Vol3{i}(:,:,2) = Wei3{i}(:,:,2) ./ ...
      (Wei3{i}(:,:,1)/pmd.Dens3{i}(1) + Wei3{i}(:,:,2)/pmd.Dens3{i}(2) + Wei3{i}(:,:,3)/pmd.Dens3{i}(3) + eps) /...
      pmd.Dens3{1}(2);
    Vol3{i}(:,:,3) = Wei3{i}(:,:,3) ./ ...
      (Wei3{i}(:,:,1)/pmd.Dens3{i}(1) + Wei3{i}(:,:,2)/pmd.Dens3{i}(2) + Wei3{i}(:,:,3)/pmd.Dens3{i}(3) + eps) /...
      pmd.Dens3{1}(3);
  end

  disp('Calculating line integrals...')
  
  if pmd.p2MD
    p2 = cell(length(Wei2), 1);
    for j = 1:length(Wei2)
      for i = 1:2
        porig2 = sinogramJ(Wei2{j}(:, :, i).*dens{pmd.tissueOrder2(j)}, degVec, r2Vec, smd.interpolation)';
        X = size(porig2, 2);
        p2{j}(:, :, i) = porig2(:,1+(X-Nr2)/2:X-(X-Nr2)/2)';
        p2{j}(:, :, i) = pixsiz * p2{j}(:, :, i);
      end
    end
  end
  
  if pmd.p3MD
    p3 = cell(length(Wei3), 1);
    for j = 1:length(Wei3)
      for i = 1:3
        porig3 = sinogramJ(Vol3{j}(:, :, i), degVec, r2Vec, smd.interpolation)';
        X = size(porig3, 2);
        p3{j}(:, :, i) = porig3(:,1+(X-Nr2)/2:X-(X-Nr2)/2)';
        p3{j}(:, :, i) = pixsiz * p3{j}(:, :, i);
      end
    end
  end
  
  % Compute monoenergetic projections
  %----------------------------------
  disp('Calculating monoenergetic projections...')
  if pmd.p2MD
    p2Low = cell(length(p2), 1);
    p2High = cell(length(p2), 1);
    for j = 1:length(p2);
      for i = 1:2
        p2Low{j}(:, :, i)  = p2{j}(:, :, i) * Cross2{pmd.tissueOrder2(j)}(1, i) * 100;
        p2High{j}(:, :, i) = p2{j}(:, :, i) * Cross2{pmd.tissueOrder2(j)}(2, i) * 100;
      end
    end
  end
  
  if pmd.p3MD
    p3Low = cell(length(p3), 1);
    p3High = cell(length(p3), 1);
    for j = 1:length(p3);
      for i = 1:3
        p3Low{j}(:, :, i)  = p3{j}(:, :, i) * pmd.Att3{pmd.tissueOrder3(j)}(1, i) * 100;
        p3High{j}(:, :, i) = p3{j}(:, :, i) * pmd.Att3{pmd.tissueOrder3(j)}(2, i) * 100;
      end
    end
  end
  
  if pmd.p2MD
    MLow = sum(sum(cat(4, p2Low{:}), 4), 3);
    MHigh = sum(sum(cat(4, p2High{:}), 4), 3);
  else
    MLow = 0;
    MHigh = 0;
  end
  
  if pmd.p3MD
    MLow = MLow + sum(sum(cat(4, p3Low{:}), 4), 3);
    MHigh = MHigh + sum(sum(cat(4, p3High{:}), 4), 3);
  end

  % Compute polychromatic projections
  % ---------------------------------
  disp('Calculating polychromatic projections...')
  if pmd.p2MD
    for i = 1:length(p2);
      if i == 1
        p = p2{i};
      else
        p = cat(3, p, p2{i});
      end
    end
  end
  if pmd.p3MD
    for i = 1:length(p3);
      if ~exist('p','var')
        p = p3{i};
      else
        p = cat(3, p, p3{i});
      end
    end
  end

  switch (useCode) % 0 = Matlab, 1 = C, 2 = OpenMP, 3 = OpenCL
    case {0, 1, 2}
      ApLow = computePolyProj(smd.ELow, uLow, smd.NLow, p, pmd.muLow);
      ApHigh = computePolyProj(smd.EHigh, uHigh, smd.NHigh, p, pmd.muHigh);
    case 3
      [ApLow, ApHigh] = computePolyProjc_opencl(smd.ELow, smd.EHigh, uLow, uHigh,...
        smd.NLow, smd.NHigh, p, pmd.muLow, pmd.muHigh);
  end
  clear('p');
  
  ZLow  = pmd.projLow + (MLow - ApLow);
  ZHigh = pmd.projHigh + (MHigh - ApHigh);
  
  % Reconstruction
  % --------------
  disp('Reconstruction...')
  recLow = reconstructIteratedProjections(ZLow, r2Vec, degVec, smd.N1, smd.dt1);
  recHigh = reconstructIteratedProjections(ZHigh, r2Vec, degVec, smd.N1, smd.dt1);
  
  pmd.recLowSet{pmd.curIterIndex} = recLow;
  pmd.recHighSet{pmd.curIterIndex} = recHigh;
  
  % Plot reconstructed maps of linear attenuation coefficients
  if iter == numbIter && gDiraPlotFigures == 1
    pmd.PlotRecLacImages(iter);
    drawnow();
  end
  
  % Tissue segmentation
  % -------------------
  
  disp('Classifying tissues...')
  [tissue2 tissue3] = tissueClassification(iter, smd, pmd);
  pmd.tissue2Set{pmd.curIterIndex} = tissue2;
  pmd.tissue3Set{pmd.curIterIndex} = tissue3;

  % Tissue decomposition
  % ---------------------------
  disp('Decomposing tissues...')
  AttE1mat = 0.01*recLow;		% Change 1/m to 1/cm
  AttE2mat = 0.01*recHigh;
  
  if pmd.p2MD
    dens = cell(length(tissue2), 1);
    Wei2 = cell(length(tissue2), 1);
    for i = 1:length(tissue2)
      [Wei2{i}, dens{i}] = MD2(AttE1mat, AttE2mat, pmd.Att2{pmd.tissueOrder2(i)},...
        pmd.Dens2{pmd.tissueOrder2(i)}, tissue2{i});
    end
  end
  
  if pmd.p3MD
    Wei3 = cell(length(tissue3), 1);
    for i = 1:length(tissue3)
      Wei3{i} = MD3(AttE1mat, AttE2mat, pmd.Att3{pmd.tissueOrder3(i)}, pmd.Dens3{pmd.tissueOrder3(i)}, tissue3{i}, 1);
    end
  end
  
  if pmd.p2MD
    pmd.densSet{pmd.curIterIndex} = dens;
    pmd.Wei2Set{pmd.curIterIndex} = Wei2;
  end
  if pmd.p3MD
    pmd.Wei3Set{pmd.curIterIndex} = Wei3;
  end

  % Plot computed mass fractions from MD2 and MD3
  if iter == numbIter && gDiraPlotFigures == 1
    if pmd.p2MD
      pmd.PlotMassFractionsFromMd2(iter);
    end
    if pmd.p3MD
      pmd.PlotMassFractionsFromMd3(iter);
    end
    drawnow();
  end
end

%% Add material decomposition of prostate to resulting data.
%
pmd.Wei3SA{1} = MD3(AttE1mat, AttE2mat, pmd.Att3SA, pmd.Dens3SA, pmd.maskSA, 1);
% Plot computed mass fractions from MD3
if gDiraPlotFigures == 1
  plotWei3(pmd.Wei3SA, pmd.name3SA);
  drawnow();
end
pmd.WeiAv = MD3SP(mean(AttE1mat(pmd.maskSA)), mean(AttE2mat(pmd.maskSA)), pmd.Att3SA, pmd.Dens3SA);
fprintf('Average mass fraction m1 = %f, m2 = %f and m3 = %f\n', pmd.WeiAv);

pmd.curIterIndex = -1; % This state variable indicates the end of DIRA.

%% Save results
save('pmd.mat', 'pmd');
save('smd.mat', 'smd');

toc
fprintf('\nDone!\n')
