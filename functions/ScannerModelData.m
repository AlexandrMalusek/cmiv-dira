classdef ScannerModelData < handle
  properties
    eL            % low x-ray tube voltage in kV
    eH            % high x-ray tube voltage in kV
    ELow          % spectrum energies for Ul
    NLow          % relative number of photons for Ul
    EHigh         % spectrum energies for Uh
    NHigh         % relative number of photons for Uh
    L             % distance source - rot. center in m
    N1            % number of detector elements after rebinning
    dt1           % detector element size = pixel distance
    interpolation % Joseph projection generation
    % Rebinning data
    N0            % number of detector elements before rebinning
    M0            % number of projections
    fact          % factor for no rebinned projections
    dfi0          % angle increment in degrees
    dfi1          % new angle increment in degrees
    M1            % new nr of projections
    dt0           % Detector element arc length [rad]
    gamma         % first angle after rebinning [rad]
    mask          % reconstruction mask
  end

  methods
    function PlotSpectra(smd, varargin)
      % Plot energy spectra
      plot(smd.EHigh, smd.NHigh, '.-', smd.ELow, smd.NLow, '.-', varargin{:})
      title('Energy spectra')
      xlabel('energy (keV)')
      ylabel('number of photons per bin')
    end
  end % methods
end % classdef
