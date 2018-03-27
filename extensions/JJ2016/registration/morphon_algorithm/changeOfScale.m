function [resDispAcc, resCertAcc] = changeOfScale(dispAcc, certAcc, ...
                                               smallSize, nextSize)
% Funktion that upsamples the accumulative displacement field and the
% accumulated certainty to the next coarsest scale. This function is a part 
% of the implementation of the Morphon.
%
% Input:    - dispAcc:      The accumulated displacement field.
%
%           - certAcc:      The accumulated certainty.
%
%           - smallSize:    The size of dispAcc and certAcc on the coarser
%                           scale.
%
%           - nextSize:     The new size of dispAcc and certAcc. 
%
%
%
% Output:   - resDispAcc:   Upscaled accumulated displacement field.
%
%           - resCertAcc:   The accumulated certainty.
%
%  See also Morphon.

resDispAcc = zeros([nextSize 3]);
resCertAcc = zeros([nextSize 3]);

scale = (nextSize./smallSize);
for dim = 1:3
    resDispAcc(:,:,:,dim) = scale(dim)*rescale3(dispAcc(:,:,:,dim), ...
                            nextSize,'linear');
                        
    resCertAcc(:,:,:,dim) = rescale3(certAcc(:,:,:,dim),nextSize,'linear');
end


end