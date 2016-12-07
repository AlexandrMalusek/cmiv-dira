function compareVolumes(firstVol,secondVol,doSuperPos,alpha,figureNum,imPlane)
% Function for comparing two slice images of two different volumes. This
% can be done in the xy, xz and yz plane. Additional the input argument
% 'doSuperPos' can be set to 1 in order to get a third image where the
% images are plotet on top of each other.
%
% input: - firstVol
%        - secondVol
%        - doSuperPos:    If the volumes should be superpositioned or not
%                         If superpos = 0 they aren't shown, if 1 they are
%        - alpha:         alpha value for the superposition
%        - imPlane:       Defines the direction of the visualization. The
%                         image planes can be shown in the 'xy', 'xz' and
%                         the 'yz' plane. 

nrImagesUsed = 2 + doSuperPos;

vol = alpha*firstVol + (1-alpha)*secondVol;



%% Show images in the XY plane
if(strcmp(imPlane, 'all') || strcmp(imPlane, 'xy'))
    % show figures
    figure(figureNum)
    subplot(1,nrImagesUsed,1)
    imagesc(firstVol(:,:,ceil(size(firstVol,3)/2)))
    colormap gray; axis 'image'
    title('original - xy plane')
    subplot(1,nrImagesUsed,2)
    imagesc(secondVol(:,:,ceil(size(secondVol,3)/2)))
    colormap gray; axis 'image'
    title('deformed - xy plane')
    
    if(doSuperPos)
        subplot(1,nrImagesUsed,3)
        imagesc(squeeze(vol(:,:,ceil(size(secondVol,3)/2))))
        colormap gray; axis 'image'
        title('Superposition - xy plane')
    end
    drawnow
end
%% Show images in the YZ plane
if(strcmp(imPlane, 'all') || strcmp(imPlane, 'yz'))
    % flip volume
    volume_yz  = permute(firstVol,[3 1 2]);
    deformed_vol_yz = permute(secondVol,[3 1 2]);
    vol_yz = permute(vol,[3 1 2]);
    
    % show figures
    figure(figureNum + 1)
    subplot(nrImagesUsed,1,1)
    imagesc(volume_yz(:,:,ceil(size(volume_yz,3)/2)))
    colormap gray; axis 'image'
    title('original - yz plane')
    subplot(nrImagesUsed,1,2)
    imagesc(deformed_vol_yz(:,:,ceil(size(deformed_vol_yz,3)/2)))
    colormap gray; axis 'image'
    title('deformed - yz plane')
    
    if(doSuperPos)
        subplot(nrImagesUsed,1,3)
        imagesc(squeeze(vol_yz(:,:,ceil(size(vol_yz,3)/2))))
        colormap gray; axis 'image'
        title('Superposition - yz plane')
    end
    drawnow
end
%% Show images in the XZ plane
if(strcmp(imPlane, 'all') || strcmp(imPlane, 'xz'))
    % flip volume
    volume_xz  = permute(firstVol,[3 2 1]);
    deformed_vol_xz = permute(secondVol,[3 2 1]);
    vol_xz = permute(vol,[3 2 1]);
    
    % show figures
    figure(figureNum + 2)
    subplot(nrImagesUsed,1,1)
    imagesc(volume_xz(:,:,ceil(size(volume_xz,3)/2)))
    colormap gray; axis 'image'
    title('original - xz plane')
    subplot(nrImagesUsed,1,2)
    imagesc(deformed_vol_xz(:,:,ceil(size(deformed_vol_xz,3)/2)))
    colormap gray; axis 'image'
    title('deformed - xz plane')
    
    if(doSuperPos)
        subplot(nrImagesUsed,1,3)
        imagesc(squeeze(vol_xz(:,:,ceil(size(vol_xz,3)/2))))
        colormap gray; axis 'image'
        title('Superposition - xz plane')
    end
    drawnow
end

end