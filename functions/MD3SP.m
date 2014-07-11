%% Three-material decomposition in DECT.
%
function [w] = MD3SP(uxl, uxh, att, dens)

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