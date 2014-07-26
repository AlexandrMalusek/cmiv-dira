classdef ScannerModelData < handle
  properties
    eL            % low x-ray tube voltage in kV
    eH            % high x-ray tube voltage in kV
    eEL           % low effective energy in keV
    eEH           % high effective energy in keV
    ELow          % spectrum energies for Ul
    NLow          % relative number of photons for Ul
    EHigh         % spectrum energies for Uh
    NHigh         % relative number of photons for Uh
    L             % distance source - rot. center in m
    alpha         % fanbeam angle in deg
    N1            % number of detector elements after rebinning
    dt1           % detector element size = pixel distance
    interpolation % Joseph projection generation
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
