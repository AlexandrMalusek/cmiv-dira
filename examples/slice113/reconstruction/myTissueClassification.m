%% Tissue classification
%
function [tissue2, tissue3] = myTissueClassification(~,~)
  %
  % Input:
  % lac1:     matrix of LACs at effective energy E1
  % lac2:     matrix of LACs at effective energy E2
  %
  % Output:
  % tissue2:  masks defining tissues decomposed using MD2
  % tissue3:  masks defining tissues decomposed using MD2
 
  % Use tissue masks derived from the original phantom.
  temp = load('tissues.mat');
  tissue2 = temp.tissue2;
  tissue3{1} = temp.tissue3{1} + temp.tissue3{2};
end
