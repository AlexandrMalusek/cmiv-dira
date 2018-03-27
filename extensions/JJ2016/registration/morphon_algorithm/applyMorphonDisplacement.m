function deformedVolume = applyMorphonDisplacement(volume,displace,interpMetod)
% Function for applying the displacement field calculated by the Morphon 
% algorithm to a volume. 
%
% Input:    - volume          The volume that we aim to deform.
%
%           - displace        The displacement field given by the 'morphon' 
%                             function. The first three dimensions are used
%                             to store the displacement fields in y, x and 
%                             z direction. The fourth dimension is used to
%                             access the y, x and z components of the 
%                             displacement field. 
%
%           - interpMetod     The method used for interpolation. Can be set
%                             to 'nearest', 'linear' and 'cubic'.
%
%
% Output:   - deformedVolume  The deformed volume. 
%
%   See also morphon.

if nargin == 2
    interpMetod = 'linear';
end

volSize = size(volume);
[y,x,z] = ndgrid(1:volSize(1),1:volSize(2), 1:volSize(3)); 

deformedVolume = ba_interp3(x,y,z,double(volume), ...
                            x - double(displace(:,:,:,2)),   ...
                            y - double(displace(:,:,:,1)),   ...
                            z - double(displace(:,:,:,3)),interpMetod);  
                                    
                                    
                                    
end