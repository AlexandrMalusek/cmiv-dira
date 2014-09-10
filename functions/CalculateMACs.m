%% Calculate mass attenuation coefficients
%
function [uMass] = CalculateMACs(moleculeStr, eE)

  lngEE = size(eE, 2);
  uMass = zeros(lngEE, 1);
  
  fId1 = fopen('data_Ar.txt', 'r');
  mass = fscanf(fId1, '%s %i %f', [4 103])';
  fclose(fId1);
  
  fId2 = fopen('data_macTable.txt', 'r');
  macs = fscanf(fId2, '%f', [101 399])';
  fclose(fId2);
  
  % String formating
    
  alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz';
  numbers = '0123456789';
  i = 1;
  while i < length(moleculeStr)
    if ismember(moleculeStr(i), alphabet) && ismember(moleculeStr(i+1), numbers)
      moleculeStr = [moleculeStr(1:i) ' ' moleculeStr(i+1:end)];
    end
    if ismember(moleculeStr(i+1), alphabet) && ismember(moleculeStr(i), numbers)
      moleculeStr = [moleculeStr(1:i) ' ' moleculeStr(i+1:end)];
    end
    if ismember(moleculeStr(i), ' ') && ismember(moleculeStr(i+1), alphabet) && ismember(moleculeStr(i+2), ' ')
      moleculeStr = [moleculeStr(1:i) '-' moleculeStr(i+1:end)];
    end
    i = i+1;
  end
  i = 1;
  if ismember(moleculeStr(1), alphabet) && ismember(moleculeStr(2), ' ')
    moleculeStr = ['-' moleculeStr(:)'];
  end
  while i < length(moleculeStr)
    if ismember(moleculeStr(i), ' ') && ismember(moleculeStr(i+1), alphabet) && ismember(moleculeStr(i+2), ' ')
      moleculeStr = [moleculeStr(1:i) '-' moleculeStr(i+1:end)];
    end
    i = i+1;
  end
    
  molecule = sscanf(moleculeStr, '%s %f', [3 inf])';
    
  lngMol = size(molecule, 1);

  atoms = zeros(lngMol, 3);
    
  for i = 1:lngMol
    atoms(i, 1:2) = mass(ismember(mass(:,1:2), molecule(i, 1:2),'rows'), 3:4);
    atoms(i, 3) = molecule(i, 3);
  end
    
  massTot = sum(atoms(:,3).*atoms(:,2));
  
  for j = 1:lngEE
    
    macsEE = macs(macs(:,1) == eE(j),2:end);
   
    uMass(j) = 0;
    for i = 1:lngMol
      uMass(j) = uMass(j) + atoms(i, 3)*(atoms(i, 2)/massTot)*macsEE(atoms(i, 1));
    end  
  end
end
