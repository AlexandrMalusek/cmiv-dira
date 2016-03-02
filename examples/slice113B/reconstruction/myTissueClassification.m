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
  airthresh  = 19;                 % Slightly below lipid   = 20 [1/m] at 50kV
  bonethresh = 33;                 % Slightly above protein = 28 [1/m]
  se = [1 1 1; 1 1 1; 1 1 1];

  % Make masks
  % ==========
  boneMask = (pmd.recLowSet{iter+1} > bonethresh);             % threshold bone
  boneMask = imdilate(boneMask,se);                            % dilate bone mask
  boneMask = 1-bwareaopen(1-boneMask, 16000);                  % fill small holes
  airMask  = smd.mask .* (pmd.recLowSet{iter+1} < airthresh);  % air mask
  ba = boneMask .* airMask;
  airMask = airMask - ba;
  softMask = smd.mask - boneMask - airMask;                    % soft tissue mask

  tissue2{1} = logical(airMask);
  tissue2{2} = logical(boneMask);
  tissue3{1} = logical(softMask);
end
