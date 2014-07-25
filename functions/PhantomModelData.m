classdef PhantomModelData < handle
  properties
    projLow       % projection for Ul
    projHigh      % projection for Uh
    projLowBH     % WBHC projection for Ul
    projHighBH    % WBHC projection for Uh
    Dens2
    Att2
    Dens3
    Att3
    muLow
    muHigh
    tissueOrder2
    tissueOrder3
    recLowSet
    recHighSet
    densSet
    Wei2Set
    Wei3Set
  end

end % classdef
