function atlasVolume = readVisualHuman
% Function for getting the VisHum atlas. The functions reads the data in
% segm_vishum in order to create an atlas volume. 
% 
% Example: 
%   atlasVolume = readVisualHuman;


fileID = fopen('segm_vishum');
header = fread(fileID, 4096,'uchar');
atlasVector = fread(fileID, 512*512*250,'uchar');
fclose(fileID);
atlasVolume = reshape(atlasVector, [512,512,250]);
atlasVolume = permute(atlasVolume(:,end:-1:1,end:-1:1),[2,1,3]);

atlasVolume = single(atlasVolume);

end