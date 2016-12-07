function [bone, rectum, prostate] = backSamlingResults(bone,...
                                                     rectum, prostate,...
                                                     originalSize)
% Function for resampling the results from the atlas segmentation in order 
% for them to have the same voxel size as the original CT volume. 
%
% Input:   - bone:         The volume containing the bones.
%
%          - rectum:       The volume containing the rectum.
%
%          - prostate:     The volume containing the prostate.
%
%          - originalSize: The size of the CT volume before the resampling.
%
%
%
% Outputs: - bone:         Resampled bone volume.
%
%          - rectum:       Resampled rectum volume.
%
%          - prostate:     Resampled prostate volume.
                                                 
[bone,    ~]   = backSampling(double(bone),     originalSize, 'cubic');
[rectum,  ~]   = backSampling(double(rectum),   originalSize, 'cubic');
[prostate,~]   = backSampling(double(prostate), originalSize, 'cubic');
end