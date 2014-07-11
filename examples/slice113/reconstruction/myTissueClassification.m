%% Tissue classification
%
function [tissue2, tissue3] = myTissueClassification(~,~)
  
  % Use tissue masks derived from the original phantom.
  temp = load('tissues.mat');
  tissue2 = temp.tissue2;
  tissue3{1} = temp.tissue3{1} + temp.tissue3{2};
end