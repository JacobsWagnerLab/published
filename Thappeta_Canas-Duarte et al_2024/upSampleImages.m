%% Resize all tif images in a selected folder by a constant scale
% -Input-
%   scale: how many times the images should be scaled
% -Output-
%   A new folder 'scaled' will be created in the selected folder, and it
%   will contain all the scaled images
function upSampleImages(scale)
srcFolder = uigetdir();
destiFolder = strcat(srcFolder, '\scaled\');

flst = dir(strcat(srcFolder, '\*.tif'));

if ~isempty(flst)
    mkdir(destiFolder);
else
    error('No tif images found in selected folder: %s', srcFolder);
end

for frame = 1:length(flst)
    im = imread(fullfile(flst(frame).folder, flst(frame).name));
    scaledImage = imresize(im, scale);
    imwrite(scaledImage, fullfile(destiFolder, strcat('scaled', flst(frame).name)));
end

end