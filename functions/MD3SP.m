%% Three-material decomposition in DECT.
%
function [w] = MD3SP(uxl, uxh, att, dens)
  %
  % Non-matrix version of MD3. Parameters uxl and uxh are scalars.
  %
  % Input:
  % uxl:  measured LAC at effective energy E1 (scalar)
  % uxh:  measured LAC at effective energy E2 (scalar)
  % att:  tabulated LACS of triplet base materials at E1 and E2
  % dens: mass density of triplet base materials
  %
  % Output:
  % w:    mass fractions of triplet base materials
  %
  % Examples:
  %
  % In the function call:
  % WeiAv = MD3SP(mean(AttE1mat(mask)), mean(AttE2mat(mask)), Att3SA, Dens3SA);
  % the following variables are used:
  %
  % >> mean(AttE1mat(mask))
  % ans =
  %     0.2720
  %
  % >> mean(AttE2mat(mask))
  % ans =
  %     0.1903
  %
  % >> Att3SA
  % Att3SA =
  %     0.2302    0.2270    1.5808
  %     0.1808    0.1772    0.4778
  %
  % >> Dens3SA
  % Dens3SA =
  %     1.0300    1.0000    1.5500
  %
  % >> WeiAv
  % WeiAv =
  %     1.0425   -0.0880    0.0455

  um1l = att(1, 1);
  um2l = att(1, 2);
  um3l = att(1, 3);
  um1h = att(2, 1);
  um2h = att(2, 2);
  um3h = att(2, 3);
  dm1 = dens(1);
  dm2 = dens(2);
  dm3 = dens(3);
  
  A = [((uxl-um1l)/dm1)-((uxl-um3l)/dm3),...
    ((uxl-um2l)/dm2)-((uxl-um3l)/dm3);...
    ((uxh-um1h)/dm1)-((uxh-um3h)/dm3),...
    ((uxh-um2h)/dm2)-((uxh-um3h)/dm3)];
  B = -[(uxl-um3l)/dm3; (uxh-um3h)/dm3];

  W = linsolve(A,B);
  w = [W(1), W(2), 1-W(1)-W(2)];
end
