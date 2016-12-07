function [atlasWithoutTableArms, bone, atlasMask] = getTissuesAtlas(atlas)
% Function that uses the labling of the VisHum atlas to remove the CT table
% and the arms from the atlas and that creates a binary volume containing  
% the bones in the atlas. 
%
% Input:    - atlas:                    The VisHum atlas
% 
%
% Output:   - atlasWithoutTableArms:    The Vishum atlas but without the
%                                       regions corresponding to the CT
%                                       table and the arms.
%
%           - bone:                     Binary volume containing the bones
%                                       in the VisHum atlas.
%
%           - atlasMask:                Binary mask that is one at the
%                                       voxels corresponding to the body. 

% Take out all bone
bone = atlas >57 & atlas <118;

% remove table
tableMask = atlas ~= 137;

% Remove arms
armMask = ~ismember(atlas, [3:4,42:43,51:56,120:121]);

atlasWithoutTableArms = armMask.*tableMask.*atlas;

% Mask showing body
atlasMask = atlasWithoutTableArms > 0;

end