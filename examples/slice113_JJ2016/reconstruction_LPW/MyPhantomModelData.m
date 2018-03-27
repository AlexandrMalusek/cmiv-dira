% Additional variables are needed in the PhantomModelData class to
% store atlas and reference images. Also, store segmentation results
% for each iteration.

classdef MyPhantomModelData < PhantomModelData
  properties
    atlasImage      % atlas image for atlas-based segmentation
    referenceImage  % reference image for histogram matching
    referenceMask   % mask containing the body in the reference image
    segments        % {1xNi cell}[Nr x Nr double] segmented tissues
  end
end % classdef
