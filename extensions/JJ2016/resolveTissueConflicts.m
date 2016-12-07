function [bone,adipose, ...
         rectum,prostate,air] = resolveTissueConflicts(bone, adipose,   ...
                                                       rectum, prostate,...
                                                       air, priority)
                                          
% Function for resolving tissue conflicts in the segmentation results. A 
% tissue conflict occurs when a voxel has been classified as more than one 
% tissue type. This function resolves the tissue conflicts by setting the 
% voxel to the tissue type with highest "priority". The priority is 
% defined by the order of the tissues in the cell array 'priority'.
%
%
% Input:    adipose:    Binary volume contaning the adipose tissue
%
%           bone:       Binary volume contaning the bone tissue
%
%           prostate:   Binary volume contaning the prostate 
%
%           rectum:     Binary volume contaning the rectum
%
%           priority:   Cell array containing the names of the tissues. 
%                       The first cell has the highest priority and the
%                       last the lowest. 
% 
% Example:
%
%  priority = {'bone', 'air' , 'prostate', 'rectum', 'adipose'};
%
%  [adipose,bone, ...
%   prostate,rectum] = resolveTissueConflicts(adipose, bone, prostate, ...
%                                             rectum, priority); 


% Order the tissues in a cell array
for k = 1:5
    switch priority{k}
        case 'bone'
             tissues{k} = bone; 
        case 'adipose'
             tissues{k} = adipose;
        case 'prostate'
             tissues{k} = prostate;
        case 'rectum'
             tissues{k} = rectum;
        case 'air'
             tissues{k} = air;
        otherwise
            error('Tissue not defined');
    end
end

% Set voxels to zero
for k = 5:-1:1
    for l = k-1:-1:1
        tissues{k}(tissues{l}) = 0;
    end
end

% Save results in original variables
for k = 1:5
    switch priority{k}
        case 'bone'
             bone = tissues{k}; 
        case 'adipose'
             adipose = tissues{k};
        case 'prostate'
             prostate = tissues{k};
        case 'rectum'
             rectum = tissues{k};
        case 'air'
             air = tissues{k};
        otherwise
            error('Tissue not defined');
    end
end



end