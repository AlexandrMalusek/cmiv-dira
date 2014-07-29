%% Tissue classification
%
function [tissue2, tissue3] = myTissueClassification(iter, smd, pmd)
  % myTissueClassification
  %
  % Input:
  % iter:    iteration number
  % smd:     scanner model data 
  % pmd:     phantom model data
  %
  % Output:
  % tissue2:  masks defining tissues decomposed using MD2
  % tissue3:  masks defining tissues decomposed using MD2
 
  % Threshold and structure element for bone tissue
  % -----------------------------------------------
  bonethresh = 33;                 % Slightly above protein = 28 [1/m]
  se = [1 1 1; 1 1 1; 1 1 1];

  % Make masks
  % ==========
  boneMask = (pmd.recLowSet{iter+1} > bonethresh); % threshold bone
  boneMask = imdilate(boneMask,se);                % dilate bone mask
  boneMask = 1-bwareaopen(1-boneMask, 4000);       % fill small holes
  softMask = smd.mask - boneMask;                  % soft tissue
  
  tissue2{1} = boneMask;
  tissue3{1} = softMask;
end
