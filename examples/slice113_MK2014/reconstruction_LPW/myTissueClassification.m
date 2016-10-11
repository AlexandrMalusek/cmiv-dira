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

  image = pmd.recLowSet{iter+1};
  [bones, adipose, prostate, muscles] =...
    segmentation(image, pmd.atlasImage, pmd.referenceImage, 10.0);

  tissue2{1} = bones;
  tissue3{1} = adipose | prostate | muscles;
  tissue3{1} = imfill(tissue3{1}, 'holes') - tissue2{1};

  pmd.segments{iter+1} = bones + 2*adipose + 4*prostate + 16*muscles;
end
