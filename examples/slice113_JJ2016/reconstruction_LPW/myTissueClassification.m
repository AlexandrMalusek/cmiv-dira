%% Tissue classification
%
function [tissue2, tissue3] = myTissueClassification(iter, smd, pmd)
  %
  % Input:
  % iter:    iteration number
  % smd:     scanner model data 
  % pmd:     phantom model data
  %
  % Output:
  % tissue2:  masks defining tissues decomposed using MD2
  % tissue3:  masks defining tissues decomposed using MD2
  
  
  
  % Set parameters for the JJ2016
  %------------------------------

  atlasVoxelSize  = [1 1 1]; % [0.91,   0.94,   5.0] 
  scaleInfo       = [1 1 1];
  nLastIterations = 2;       % The registration is performed for the last 
                             % nLastIterations iterations only to save time.
  
  % Start at iteration
  nrIterations = max(pmd.savedIter(:));
  performRegistration = iter + nLastIterations - nrIterations > 0;
  

  % Extend volumes to 3D
  %---------------------
  sizeInZ = 50;     % The number of slices in z direction.
  
  volume        = repmat( pmd.recLowSet{iter+1},[1,1,sizeInZ]);
  referenceVol  = repmat( pmd.referenceImage,   [1,1,sizeInZ]);
  referenceMask = repmat( pmd.referenceMask,    [1,1,sizeInZ]);
  atlas         = repmat( pmd.atlasImage,       [1,1,sizeInZ]);
  
  % Segmentation
  %-------------
  
 [compactBone, ...
  bone,      ...
  adipose, rectum,   ...
  prostate, air,     ...
  remainingTissues]   = segmentation(volume,scaleInfo,referenceVol, ...
                                     referenceMask, atlas, ...
                                     atlasVoxelSize, performRegistration);
  
  % Select the central slice of the volume
  %---------------------------------------
  centerSlice = ceil(sizeInZ/2);
  
  compactBone       = compactBone(:,:,centerSlice); 
  bone              = bone(:,:,centerSlice); 
  adipose           = adipose(:,:,centerSlice); 
  rectum            = rectum(:,:,centerSlice); 
  prostate          = prostate(:,:,centerSlice); 
  air               = air(:,:,centerSlice); 
  remainingTissues  = remainingTissues(:,:,centerSlice); 
 
  
  % Save results for DIRA
  %-----------------------
  tissue2{1} = bone;
  tissue3{1} = adipose | prostate | rectum | remainingTissues;

  pmd.segments{pmd.curIterIndex} = uint16(bone + 2*adipose + 4*rectum  ...
                                        +  8*prostate + 16*air         ...
                                        + 32*remainingTissues); 
                               
end
