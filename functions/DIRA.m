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

%% Load data and initialize variables
% ------------------------------------
disp('Loading data and initializing variables...')

% Projection data.
sinograms = load(sinogramsFileName);
projLow = sinograms.projLow;
projHigh = sinograms.projHigh;

% Projection data with water beam hardening correction
sinogramsBH = load(sinogramsBhFileName);
projLowBH = sinogramsBH.projLowBH;
projHighBH = sinogramsBH.projHighBH;

% Spectrum
spectra =  load(spectraFileName);
ELow      = spectra.currSpectLow(1:75, 1);  % spectrum energies
NLow      = spectra.currSpectLow(1:75, 2);  % relative number of photons
EHigh   = spectra.currSpectHigh(1:135, 1);  % spectrum energies
NHigh   = spectra.currSpectHigh(1:135, 2);  % relative number of photons

% Prostate mask
load(prostateMaskFileName);


%% CT scan geometry and Joseph metod
% ------------------------------------
gLen = alpha * pi/ 180;          % detector arc length [rad]

% Rebinning of fan beam data
%---------------------------
gamma = atan(0.5*gLen/L);        % first angle after rebinning [rad]
firstang = gamma;                % first angle after rebinning [rad]

% Init data for Joseph projection generation
% ------------------------------------------
Nr = N1;                         % number of detector elements
interpolation = 2;
Nphi = size(projLow, 2);         % number of projection angles after rebinning
degLen = 180; 		         % angular interval in degrees
pixsiz = dt1;                    % pixel size
Nr2 = 2*floor((Nr-1)/2)+1;       % resolution of detector elements
degVec = ((0:-1:(-Nphi+1))*(degLen/Nphi))- gamma * 180/pi;
r2Vec  = (-(Nr2-1)/2:1:(Nr2-1)/2);

% Compute I0 for both spectra
% --------------------------
sp = zeros(1,78);
for k=2:length(ELow)-1;                           
  sp(k) = (ELow(k)*NLow(k)*(ELow(k+1)-ELow(k-1)));
end
uLow = sum(sp)/2;
sp = zeros(1,138);
for k=2:length(EHigh)-1;                           
  sp(k) = (EHigh(k)*NHigh(k)*(EHigh(k+1)-EHigh(k-1)));
end
uHigh = sum(sp)/2;

%% Reconstruction No.0
%
fprintf('\nStarting initial reconstruction...\n')

phm1 = iradon(projLowBH, degVec, 'linear', 'Hann', 1, N1)/dt1;
phm2 = iradon(projHighBH, degVec, 'linear', 'Hann', 1, N1)/dt1;

recLowSet = cell(numbiter+1, 1);
recHighSet = cell(numbiter+1, 1);
recLowSet{1} = phm1;
recHighSet{1} = phm2;

plotPolyImgs(phm1, phm2, 0)
drawnow();

%% Inital Tissue segmentation
% 
disp('Classifying tissues...')
[tissue2 tissue3] = feval(tissueClass, phm1, phm1);

% Tissue decomposition
disp('Decomposing tissues...')
AttE1mat = 0.01*phm1;		   % Change 1/m to 1/cm
AttE2mat = 0.01*phm2;

if p2MD
  dens = cell(length(tissue2), 1);
  Wei2 = cell(length(tissue2), 1);
  for i = 1:length(tissue2)
    [Wei2{i}, dens{i}] = MD2(AttE1mat, AttE2mat, Att2{tissueOrder2(i)},...
      Dens2{tissueOrder2(i)}, tissue2{i});
  end
end

if p3MD
  Wei3 = cell(length(tissue3), 1);
  for i = 1:length(tissue3)
    Wei3{i} = MD3(AttE1mat, AttE2mat, Att3{tissueOrder3(i)},...
      Dens3{tissueOrder3(i)}, tissue3{i});
  end
end

densSet = cell(numbiter+1, 1);
Wei2Set = cell(numbiter+1, 1);
Wei3Set = cell(numbiter+1, 1);
if p2MD
  densSet{1} = dens;
  Wei2Set{1} = Wei2;
end
if p3MD
  Wei3Set{1} = Wei3;
end

if p2MD
  plotWei2Dens(Wei2, dens)
end
if p3MD
  plotWei3(Wei3);
end

drawnow();

%% Iterate
%
iterno = numbiter;
for iter = 1:numbiter
  % Projection generation with Joseph
  %----------------------------------
  % For every 2 or 3 MD added this part have to be extended with a loop
  % for projection calculation, monoenergetic loop and summation and
  % polychromatic concatination.
  
  fprintf('\nStarting iteration %d...\n', iter);
  
  disp('Calculating line integrals...')
  
  if p2MD
    p2 = cell(length(Wei2), 1);
    for j = 1:length(Wei2)
      for i = 1:2
        porig2 = sinogramJ(Wei2{j}(:, :, i).*dens{tissueOrder2(j)}, degVec, r2Vec, interpolation)';
        X = size(porig2, 2);
        p2{j}(:, :, i) = porig2(:,1+(X-Nr2)/2:X-(X-Nr2)/2)';
        p2{j}(:, :, i) = pixsiz * p2{j}(:, :, i);
      end
    end
  end
  
  if p3MD
    p3 = cell(length(Wei3), 1);
    for j = 1:length(Wei3)
      for i = 1:3
        porig3 = sinogramJ(Wei3{j}(:, :, i), degVec, r2Vec, interpolation)';
        X = size(porig3, 2);
        p3{j}(:, :, i) = porig3(:,1+(X-Nr2)/2:X-(X-Nr2)/2)';
        p3{j}(:, :, i) = pixsiz * p3{j}(:, :, i);
      end
    end
  end
  
  % Compute monoenergetic projections
  %----------------------------------
  disp('Calculating monoenergetic projections...')
  if p2MD
    p2Low = cell(length(p2), 1);
    p2High = cell(length(p2), 1);
    for j = 1:length(p2);
      for i = 1:2
        p2Low{j}(:, :, i)  = p2{j}(:, :, i) * Cross2{tissueOrder2(j)}(1, i) * 100;
        p2High{j}(:, :, i) = p2{j}(:, :, i) * Cross2{tissueOrder2(j)}(2, i) * 100;
      end
    end
  end
  
  if p3MD
    p3Low = cell(length(p3), 1);
    p3High = cell(length(p3), 1);
    for j = 1:length(p3);
      for i = 1:3
        p3Low{j}(:, :, i)  = p3{j}(:, :, i) * Att3{tissueOrder3(j)}(1, i) * 100;
        p3High{j}(:, :, i) = p3{j}(:, :, i) * Att3{tissueOrder3(j)}(2, i) * 100;
      end
    end
  end
  
  if p2MD
    MLow = sum(sum(cat(4, p2Low{:}), 4), 3);
    MHigh = sum(sum(cat(4, p2High{:}), 4), 3);
  else
    MLow = 0;
    MHigh = 0;
  end
  
  if p3MD
    MLow = MLow + sum(sum(cat(4, p3Low{:}), 4), 3);
    MHigh = MHigh + sum(sum(cat(4, p3High{:}), 4), 3);
  end

  % Compute polychromatic projections
  % ---------------------------------
  disp('Calculating polychromatic projections...')
  if p2MD
    for i = 1:length(p2);
      if i == 1
        p = p2{i};
      else
        p = cat(3, p, p2{i});
      end
    end
  end
  if p3MD
    for i = 1:length(p3);
      if ~exist('p','var')
        p = p3{i};
      else
        p = cat(3, p, p3{i});
      end
    end
  end

  ApLow = computePolyProj(ELow, uLow, NLow, p, muLow);
  ApHigh = computePolyProj(EHigh, uHigh, NHigh, p, muHigh);
  clear('p');
  
  ZLow  = projLow + (MLow - ApLow);
  ZHigh = projHigh + (MHigh - ApHigh);
  
  % Reconstruction
  % --------------
  disp('Reconstruction...')
  recLow = iradon(ZLow, degVec, 'linear', 'Hann', 1, N1)/dt1;
  recHigh = iradon(ZHigh, degVec, 'linear', 'Hann', 1, N1)/dt1;
  
  recLowSet{iter+1} = recLow;
  recHighSet{iter+1} = recHigh;
  
  if iter == numbiter
    plotPolyImgs(recLow, recHigh, iter);
    drawnow();
  end
  
  % Tissue segmentation
  % -------------------
  
  disp('Classifying tissues...')
  [tissue2 tissue3] = feval(tissueClass);
  
  % Tissue decomposition
  % ---------------------------
  disp('Decomposing tissues...')
  AttE1mat = 0.01*recLow;		% Change 1/m to 1/cm
  AttE2mat = 0.01*recHigh;
  
  if p2MD
    dens = cell(length(tissue2), 1);
    Wei2 = cell(length(tissue2), 1);
    for i = 1:length(tissue2)
      [Wei2{i}, dens{i}] = MD2(AttE1mat, AttE2mat, Att2{tissueOrder2(i)},...
        Dens2{tissueOrder2(i)}, tissue2{i});
    end
  end
  
  if p3MD
    Wei3 = cell(length(tissue3), 1);
    for i = 1:length(tissue3)
      Wei3{i} = MD3(AttE1mat, AttE2mat, Att3{tissueOrder3(i)}, Dens3{tissueOrder3(i)}, tissue3{i});
    end
  end
  
  if p2MD
    densSet{iter+1} = dens;
    Wei2Set{iter+1} = Wei2;
  end
  if p3MD
    Wei3Set{iter+1} = Wei3;
  end
  if iter == numbiter
    if p2MD
      plotWei2Dens(Wei2, dens)
    end
    if p3MD
      plotWei3(Wei3);
    end
  end
  drawnow();
end

%% Add material decomposition of prostate to resulting data.
%
% Define the prostate mask
mask = maskProst;

Wei3SA{1} = MD3(AttE1mat, AttE2mat, Att3SA, Dens3SA, mask);
Wei3Set{iter+2} = Wei3SA;
plotWei3(Wei3SA);
WeiAv = MD3SP(mean(AttE1mat(mask)), mean(AttE2mat(mask)), Att3SA, Dens3SA);
fprintf('Average mass fraction m1 = %f, m2 = %f and m3 = %f\n', WeiAv);

%% Save results
save(resultsFileName, 'recLowSet', 'recHighSet', 'densSet',...
  'Wei2Set', 'Wei3Set');

toc
fprintf('\nDone!\n')