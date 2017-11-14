% Ne: number of elements (= 103)
%
% Note: data_macTable.txt contains data for elements with Z = 1..100

classdef Material < handle
  properties
    nameStr         % material name
    compositionStr  % mixture defined by the number of atoms, e.g. 'H2O1' for water
    compositionType % atomic fraction ('atFr') or mass fraction ('maFr')
    density         % mass density in g/cm^3
    Ar              % [Ne x 1 double] vector of relative atomic masses
    W               % [Ne x 1 double] vector of elemental mass fractions
    chemSymbol      % [1 x Ne cell] cell array of chemical symbols
  end

  methods
  function matObj = Material(nameStr, density, compositionStr, compositionType)
    % Constructor
    %
    % nameStr         material name
    % density         mass density in g/cm^3
    % compositionStr  chemical composition, for instance 'H0.33O0.67'
    % compositionType compositionStr contain atomic fractions ('atFr') or
    %                 mass fractions ('maFr')
    %
    % If compositionType is omitted, 'atFr' is assumed for backward compatibility

    % Initialization
    if nargin == 0
      matObj.nameStr = '';
    else
      matObj.nameStr = nameStr;
    end
    if nargin <= 1
      matObj.density = 0.0;
    else
      matObj.density = density;
    end
    if nargin <= 2
      matObj.compositionStr = '';
    else
      matObj.compositionStr = compositionStr;
    end
    if nargin <= 3
      matObj.compositionType = 'atFr';
    else
      matObj.compositionType = compositionType;
    end
    matObj.Ar = zeros(103, 1);
    matObj.W = zeros(103, 1);

    % Set mass fractions W and relative atomic masses Ar
    switch matObj.compositionType
      case 'atFr'
        setMaFrFromAtFr(matObj, compositionStr);
      case 'maFr'
        setMaFrFromMaFr(matObj, compositionStr);  
    end

    matObj.chemSymbol = {...
     'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg',...
     'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V',...
     'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se',...
     'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh',...
     'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba',...
     'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho',...
     'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt',...
     'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac',...
     'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',...
     'Md', 'No', 'Lr'};
  end
  
  function setMaFrFromAtFr(matObj, atFrStr)
    % Set normalized mass fractions
    atoms = fractionsFromStr(matObj, atFrStr);
    massTot = sum(atoms(:,3));
    nAtoms = size(atoms,1);
    massTot = sum(atoms(:,3).*atoms(:,2));
    for i = 1:nAtoms
      matObj.W(atoms(i,1)) = atoms(i,2) * atoms(i,3) / massTot;
      matObj.Ar(atoms(i,1)) = atoms(i,2);
    end
  end
  
  function setMaFrFromMaFr(matObj, maFrStr)
    % Set normalized mass fractions
    atoms = fractionsFromStr(matObj, maFrStr);
    massTot = sum(atoms(:,3));
    nAtoms = size(atoms,1);
    for i = 1:nAtoms
      matObj.W(atoms(i,1)) = atoms(i,3) / massTot;
      matObj.Ar(atoms(i,1)) = atoms(i,2);
    end
  end
  
  function atoms = fractionsFromStr(matObj, compositionStr)
    % Return a matrix containing fractions
    %
    % Examples of correct compositionStr:
    % 'H2O1': integer fractions
    % 'H2.0O1.0': floating point fractions
    % 'H2.0e+01O1.0e+01': floating point fractions in scientific notation
    % 'H-2.0O-1.0': negative fractions

    fId1 = fopen('data_Ar.txt', 'r');
    tabAr = fscanf(fId1, '%s %i %f', [4 103])';
    fclose(fId1);
    
    % Replace 'E' or 'e' in numbers with '#'
    s1 = regexprep(compositionStr,'([^A-Z])([eE])([+-]*)(\d)','$1#$3$4');
    
    % Add spaces arround chemical symbols
    s1 = regexprep(s1,'([a-zA-Z])(\d)','$1 $2'); % between letter - number
    s1 = regexprep(s1,'(\d)([a-zA-Z])','$1 $2'); % between number - letter
    s1 = regexprep(s1,'([a-zA-Z])(-)','$1 $2'); % between letter - minus sign (a negative fraction)
    
    % Add '-' in front of single letter chemical symbols
    s1 = regexprep(s1,'[ ]([a-zA-Z])[ ]',' -$1 '); % between spaces
    s1 = regexprep(s1,'^([a-zA-Z])[ ]',' -$1 ');   % first letter
    
    % Replace '#' with 'E'
    s1 = regexprep(s1,'#','E');
    
    molecule = sscanf(s1, '%s %f', [3 inf])';
    lenMolecule = size(molecule, 1);
    
    atoms = zeros(lenMolecule, 3);
    
    for i = 1:lenMolecule
      atoms(i, 1:2) = tabAr(ismember(tabAr(:,1:2), molecule(i, 1:2),'rows'), 3:4);
      atoms(i, 3) = molecule(i, 3);
    end
  end
  
  function mac = computeMac(matObj, energy)
  % Compute mass attenuation coefficients (MACs) at specified energies
  %
  % energy    number or a vector with photon energies in keV
  % mac       MAC in cm^2/g
  %
  % Examples:
  % mat = Material('water', 1.0, 'H2O1');
  % mat.computeMac(30.0); mat.computeMac([20, 50]');

  persistent macTab;

  if isempty(macTab)
    macTab =  dlmread('data_macTable.txt');
  end

  % Calculate MAC vector for the material
  macVec = macTab(:,2:end) * matObj.W(1:100);

  % Calculate MAC at specified energies

  % Linear interpolation
  % mac = interp1(macTab(:,1), macVec, energy);

  % Linear interpolation in log-log coordinates
  mac = exp(interp1(log(macTab(:,1)), log(macVec), log(energy)));

  end

  function lac = computeLac(matObj, energy)
  % Compute linear attenuation coefficients (LACs) at specified energies
  %
  % energy    number or a vector with photon energies in keV
  % lac       LAC in 1/cm

  lac = matObj.density * matObj.computeMac(energy);
  end

  function meac = computeMeac(matObj, energy)
  % Compute mass energy absorption coefficients (MEACs) at specified energies
  %
  % energy    number or a vector with photon energies in keV
  % mac       MEAC in cm^2/g

  persistent meacTab;

  if isempty(meacTab)
    meacTab =  dlmread('data_meacTable.txt');
  end

  % Calculate MEAC vector for the material
  meacVec = meacTab(:,2:end) * matObj.W(1:92);

  % Calculate MEAC at specified energies
  % Linear interpolation in log-log coordinates, MeV -> keV
  meac = exp(interp1(log(1000*meacTab(:,1)), log(meacVec), log(energy)));
  end

  function elemMasFra = computeElemMasFra(matObj, Z)
    % Compute elemental mass fraction of an element with the atomic number Z 
    %
    % Z:  atomic number

    elemMasFra = matObj.W(Z);
  end

  function elemMassFra = computeElemMassFraVec(matObj)
    % Compute elemental mass fractions for the material. Sum of the fractions equals 1.

    elemMassFra = matObj.W;
  end

  function elemAtomFra = computeElemAtomFraVec(matObj)
    % Compute elemental atomic fractions for the metrial. Sum of the fractions equals 1.
    
    elemAtomFra = zeros(103,1);
    sel = matObj.W ~= 0; % indices for nonzero mass fractions
    elemAtomFra(sel) = matObj.W(sel) ./ matObj.Ar(sel);
    elemAtomFra = elemAtomFra / (sum(elemAtomFra)); % Normalize
  end

  function elemElecFra = computeElemElecFraVec(matObj)
    % Compute elemental electronic fractions for the metrial. Sum of the fractions equals 1.

    elemElecFra = zeros(103,1);
    sel = matObj.W ~= 0; % indices for nonzero mass fractions
    Z = (1:103)';  % column vector with atomic numbers
    elemElecFra(sel) = matObj.W(sel) .* Z(sel) ./ matObj.Ar(sel);
    elemElecFra = elemElecFra / (sum(elemElecFra)); % Normalize
  end

  end % methods
end % classdef
