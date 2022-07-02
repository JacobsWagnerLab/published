function [imageStack] = loadIMstack(fname)
info = imfinfo(fname);
imageStack = [];

for k = 1:length(info)
    currentImage = imread(fname, k, 'Info', info);
    imageStack(:,:,k) = currentImage;
end 
end