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
    fId1 = fopen('data_Ar.txt', 'r');
    tabAr = fscanf(fId1, '%s %i %f', [4 103])';
    fclose(fId1);
    
    % Replace 'E' or 'e' in numbers with '#'
    s1 = regexprep(compositionStr,'([^A-Z])([eE])([+-]*)(\d)','$1#$3$4');
    
    % Add spaces arround chemical symbols
    s1 = regexprep(s1,'([a-zA-Z])(\d)','$1 $2'); % between letter - number
    s1 = regexprep(s1,'(\d)([a-zA-Z])','$1 $2'); % between number - letter
    
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

  end % methods
end % classdef
