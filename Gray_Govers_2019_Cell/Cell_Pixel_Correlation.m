function [corrMatrix, mask, valsOut, valsOut2] = Cell_Pixel_Correlation(images, cellList, Frame, Cell, dx_from_center, pole_length)
%{
-About-
Performs a simple correlation (see MATLAB's help on corr/corrcoef) on pixels
from the images in 'images'. cellPxCorr identifies pixels in the biological
cell Cell and performs correlation on the pixel from all images. A
correlation matrix is returned. 

-Inputs-
images: a cell array of images

cellList: a cellList generated from Oufti that corresponds to input images

Frame: the cellList Frame that corresponds to input images

Cell: the biological cell number in the cellList to perform correlation on

dx_from_center = 1;
This determines how many pixels away from the centerline you want to
calculate the correlation. Must add up to an integer. A value of one
means that the correlation area will be two pixels wide

pole_length = 8;
This is the number of pixels away from the pole that you want to
calculate the correlation for. This is to avoid the artificial positive
correlation at the pole that results from the decrease in volume.

-Outputs-
This function performs pixel correlation on two signals of interest in a
cellList and returns the correlations (see MATLAB corr).


-Keywords-
Pixel correlation
Ribosome segregation

-Dependencies-
Extract_Cell_Pixels.m
Cell_Projection.m
Taylor_Smooth.m

-Author-
Bradley Parry October 10, 2018
%}

%Set variables
corrMatrix = [];
mask = [];
valsOut = [];
valsOut2 = [];

%If no parameters are specified, use default parameters. These parameters
%are arbitrary and it is recommended to extensively scan parameters for
%each given bacterial species and growth condition to confirm that you are
%using the best possible parameters.
if nargin == 4
    dx_from_center = 2.5;
    pole_length = 6;
end

%Find the correct cell entry based on the cellId from the Oufti cellList
ix = find(cellList.cellId{Frame} == Cell);
if isempty(ix)
    disp('cellId not found in cellList.cellId{Frame}')
    return
else
    Cell = ix;
end

%Check to make sure the cell is a real cell, not an accidental output from
%Oufti
if ~isfield(cellList.meshData{Frame}{Cell},'mesh') || length(cellList.meshData{Frame}{Cell}.mesh) < 6
    return
end

%Set additional variables
cellMesh = cellList.meshData{Frame}{Cell}.mesh;
imSize = size(images{1});

%Determine the pixels to be used for the correlation calculation
[imageInds, mask] = Extract_Cell_Pixels(cellMesh, imSize, dx_from_center, pole_length);

%Calculate the correlation using the image coordinates/pixels that were
%extracted using the Extract_Cell_Pixels function
for k = 1:length(images)

    vals = images{1}(imageInds);
    vals2 = images{2}(imageInds);
    corrMatrix = corrcoef(vals,vals2);
    
    valsOut = vals;
    valsOut2 = vals2;

end

return