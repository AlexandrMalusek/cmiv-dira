% Test selected functions in Material.m. A failed test reports 'failed',
% otherwise 'OK' is reported.
%
% Usage:
% >> t_Material
% 001: OK
% ...
% 010: OK

geps = 1e-10; % Global epsilon

% 001 Test the default calculation of mass fractions
mat = Material('water', 1.0, 'H2O1');
t_W = [0.111898344072365, 0, 0, 0, 0, 0, 0, 0.888101655927635, 0, 0]';
if (sum(abs(t_W - mat.W(1:10)) < geps) == 10)
  disp('001: OK');
else
  disp('001: failed');
end

% 002 Test the default calculation of mass fractions
mat = Material('water', 1.0, 'H2O1', 'atFr');
t_W = [0.111898344072365, 0, 0, 0, 0, 0, 0, 0.888101655927635, 0, 0]';
if (sum(abs(t_W - mat.W(1:10)) < geps) == 10)
  disp('002: OK');
else
  disp('002: failed');
end

% 003 Test the default calculation of mass fractions for negative mass fractions
mat = Material('water', 1.0, 'H-2O1', 'atFr');
t_W = [-0.144161126812133, 0, 0, 0, 0, 0, 0, 1.144161126812133, 0, 0]';
if (sum(abs(t_W - mat.W(1:10)) < geps) == 10)
  disp('003: OK');
else
  disp('003: failed');
end

% 004 Test the default calculation of mass fractions for negative mass fractions
mat = Material('water', 1.0, 'H-2O-1', 'atFr');
t_W = [0.111898344072365, 0, 0, 0, 0, 0, 0, 0.888101655927635, 0, 0]';
if (sum(abs(t_W - mat.W(1:10)) < geps) == 10)
  disp('004: OK');
else
  disp('004: failed');
end

% 005 Test the default calculation of mass fractions for exponential numbers
mat = Material('water', 1.0, 'H-2.0e-01O-1.0e-01', 'atFr');
t_W = [0.111898344072365, 0, 0, 0, 0, 0, 0, 0.888101655927635, 0, 0]';
if (sum(abs(t_W - mat.W(1:10)) < geps) == 10)
  disp('005: OK');
else
  disp('005: failed');
end

% 006 Test the calculation of mass fractions for 'maFr'
mat = Material('test', 1.0, 'Fe0.04E+3Xe0.06E+3', 'maFr');
t_W = [0.4, 0.6]';
if (sum(abs(t_W - mat.W([26, 54])) < geps) == 2)
  disp('006: OK');
else
  disp('006: failed');
end

% 007 Test computeMac(matObj, energy), scalar
mat = Material('water', 1.0, 'H2O1');
mac = mat.computeMac(70.0);
t_mac = 1.929E-01; % From NIST's XCOM
if (abs(t_mac - mac) < 0.001)
  disp('007: OK');
else
  disp('007: failed');
end

% 008 Test computeMac(matObj, energy), column vector
mat = Material('water', 1.0, 'H2O1');
mac = mat.computeMac([33.5, 88.5]');
t_mac = [3.239E-01, 1.775E-01]'; % cm^2/g, from NIST's XCOM
if (sum(abs(t_mac - mac) < 0.001) == 2)
  disp('008: OK');
else
  disp('008: failed');
end

% 009 Test computeLac(matObj, energy), column vector
mat = Material('water', 2.0, 'H2O1'); % Double the density
lac = mat.computeLac([33.5, 88.5]');
t_lac = [0.6478, 0.3550]'; % 1/cm, NIST's XCOM values multiplied by 2.0 g/cm^3
if (sum(abs(t_lac - lac) < 0.001) == 2)
  disp('009: OK');
else
  disp('009: failed');
end

% 010 Test computeMeac(matObj, energy), column vector
mat = Material('water', 1.0, 'H2O1');
meac = mat.computeMeac([30.0, 80.0]');
t_meac = [1.557E-01, 2.597E-02]'; % cm^2/g, from NIST's XCOM
if (sum(abs(t_meac - meac) < 0.001) == 2)
  disp('010: OK');
else
  disp('010: failed');
end
