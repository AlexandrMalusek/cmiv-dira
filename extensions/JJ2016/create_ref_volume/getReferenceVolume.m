function [newRef, refMask] = getReferenceVolume(ref)
% Function used for creating the reference volume for the histogram 
% matching. 
%
% Input:            - ref:      3D volume containing the CT data.
%
% Output:           - newRef:   The reference volume with changed contrast.
%
%                   - refMask:  Binary volume containing the body of the
%                               patient.

newRef = change_contrast(ref); 

% Create a mask showing the body.
refMask = newRef > 0;

% Find biggest object - This results in that the table is removed
labeled = bwlabeln(refMask,6);

tmp     = labeled(refMask > 0);
M       = hist(tmp(:),max(tmp(:)));
refMask = labeled == (find(max(M) == M));

newRef = newRef.*refMask;

end