function [imageStack] = loadIMstack(fname)
%load the images found in the file fname and compile them into a
%3dimensional matrix wherethe output is organized as
%   imageStack(imageRows,imageColumns,frames)
%
%Brad Parry, Christine Jacobs-Wagner lab; April 2016

info = imfinfo(fname);
imageStack = [];

for k = 1:length(info)
    currentImage = imread(fname, k, 'Info', info);
    imageStack(:,:,k) = currentImage;
end 
end