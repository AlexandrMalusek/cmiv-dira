function mipIm = maximumImageProjection(volume,directionOfProjection)
% Function that creates maximum intensity projection (MIP) image, where the
% values are projected along the 'xy', 'xz' or the 'yz' axis. 
%
% Input: - volume:                Volume containing gray values.
%
%        - directionOfProjection: The direction of the projection. Can be
%                                 set to 'x', 'y' or 'z'.
%
% Output: - mipIm:                The MIP image.

switch directionOfProjection
    case 'y'
        mipIm = max(volume,[],1);
        mipIm = permute(mipIm,[3 2 1]);
    case 'x'
        mipIm = max(volume,[],2);
        mipIm = permute(mipIm,[3 1 2]);
    case 'z'
        mipIm = max(volume,[],3);
    otherwise
        error('No okay direction!')
end


end