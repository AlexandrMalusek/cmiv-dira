function weight = atlasWeightFunction(atlas,slopeSize,filterDir)
% This funktion creates a weight funtion with the same size as the atlas.
% The weight consists of a central part that has the value one and two 
% slopes, one on either side of the central part, in the direction  
% specified by the filter direction 'filterDir'. 
% 
%
% Input:    atlas:         - The atlas volume.
%
%           slopeSize:     - Integer specifying the size of the two slopes.
%
%           filterDir:     - Direction of the filter. Can be set to 'x',
%                            'y' and 'z'.
%
% Output:   weight:        - A weight function, with the same size as the
%                            atlas volume. 

[sy,sx,sz] = size(atlas);

% Gaussian function
gaussFunction = fspecial('gaussian',[1 2*slopeSize],(2*slopeSize)/6);
gaussFunction = gaussFunction/max(gaussFunction);

% Logarithmic function
logFunction = log(1:ceil(slopeSize));
logFunction = logFunction/max(logFunction);

%% Create the weight
switch filterDir
    case 'z'
        filterSize = sz;
        zFilter = [gaussFunction(1:ceil(slopeSize)).*logFunction ...
                  ones(1,filterSize-2*slopeSize),  ...
                  gaussFunction((slopeSize+1):end).*logFunction(end:-1:1)];
        
        shiftedWeight = permute(zFilter,[1,3,2]);
        weight         = single(repmat(shiftedWeight,[sy,sx,1]));
    case 'x'
        filterSize = sx;
        xFilter = [gaussFunction(1:ceil(slopeSize)).*logFunction ...
                  ones(1,filterSize-2*slopeSize), ...
                  gaussFunction((slopeSize+1):end).*logFunction(end:-1:1)];
        
        weight  = single(repmat(xFilter,[sy,1,sz]));       
    case 'y' 
        filterSize = sy;
        yFilter = [gaussFunction(1:ceil(slopeSize)).*logFunction ...
                  ones(1,filterSize-2*slopeSize),   ...
                  gaussFunction((slopeSize+1):end).*logFunction(end:-1:1)];

        shiftedWeight = permute(yFilter,[2 1 3]);
        weight        = single(repmat(shiftedWeight,[1,sx,sz]));             
    otherwise
        
end

end

