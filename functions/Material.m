% Ne: number of elements (= 103)

classdef Material
  properties
    nameStr     % material name
    moleculeStr % mixture defined by the number of atoms, e.g. 'H2O1' for water 
    density     % mass density in g/cm^3
    Ar          % [Ne x 1 double] vector of relative atomic masses
    W           % [Ne x 1 double] vector of elemental mass fractions
  end
  
  methods
    function matObj = Material(nameStr, density, moleculeStr)
      % Constructor
      matObj.nameStr = nameStr;
      matObj.moleculeStr = moleculeStr;
      matObj.density = density;
      matObj.W = zeros(103, 1);

      fId1 = fopen('data_Ar.txt', 'r');
      mass = fscanf(fId1, '%s %i %f', [4 103])';
      fclose(fId1);
      matObj.Ar = mass(:, 4);

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
      
      for i = 1:lngMol
	  matObj.W(atoms(i,1)) = atoms(i,2) * atoms(i,3) / massTot;
      end
    end
    
  end % methods
end % classdef
