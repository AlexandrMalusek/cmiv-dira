function [bone,adipose, ...
         muscles,prostate] = resolveTissueConflicts(bone,adipose,...
                                                 muscles,prostate,priority)
% Function for resolving tissue conflicts in the segmentation results. A 
% tissue conflict occurs when a pixel has been classified as more than one 
% tissue type. This function resolves the tissue conflicts by setting the 
% pixel to the tissue type with highest 'priority'. The priority is 
% defined by the order of the tissues in the cell array priority.
%
% Input:    adipose:    Binary image contaning the adipose tissue
%
%           bone:       Binary image contaning the bone tissue
%
%           prostate:   Binary image contaning the prostate 
%
%           muscles:    Binary image contaning the muscles
%
%           priority:   Cell array containing the names of the tissues. 
%                       The first cell has the highest priority and the
%                       last the lowest. 
% 
% Example:
%
%  priority = {'prostate', 'bone', 'muscles', 'adipose'};
%
%  [adipose,bone,...
%  prostate,rectum] = resolveTissueConflicts(adipose,bone,prostate, ...
%                                            rectum,priority); 


% Order the tissues in a cell array
for k = 1:4
    switch priority{k}
        case 'bone'
             tissues{k} = bone; 
        case 'adipose'
             tissues{k} = adipose;
        case 'prostate'
             tissues{k} = prostate;
        case 'muscles'
             tissues{k} = muscles;
        otherwise
            error('Tissue not defined');
    end
end

% Set voxels to zero
for k = 4:-1:1
    for l = k-1:-1:1
        tissues{k}(tissues{l}) = 0;
    end
end

% Save results in original variables
for k = 1:4
    switch priority{k}
        case 'bone'
             bone = tissues{k}; 
        case 'adipose'
             adipose = tissues{k};
        case 'prostate'
             prostate = tissues{k};
        case 'muscles'
             muscles = tissues{k};
        otherwise
            error('Tissue not defined');
    end
end



end