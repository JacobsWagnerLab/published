function [imageInds, mask] = Extract_Cell_Pixels(cellMesh, imSize, dC, poleLen)
%{
-About-
This function uses Cell_Projection.m to extract pixels within dC of the 
cell centerline and greater than poleLen from the pole.

-Inputs-
cellMesh = cellList.meshData{1}{1}.mesh;

imSize = imread(files{1}{1});

poleLen = 6; %pixels within this distance of the pole will be ignored

dC = 3; %radius from the centerline that pixels are considered within

-Outputs-
imageInds:  extracted pixel indices with respect to image.

mask:  mask of the extracted indices.


-Keywords-
Pixel correlation
Ribosome segregation

-Dependencies-
Cell_Projection.m
Taylor_Smooth.m

-Author-
Bradley Parry October 10, 2018
%}

%convert a cell mesh into a polygon outline
poly = double(cat(1,cellMesh(:,1:2),flipud(cellMesh(:,3:4)))); 
%use the cell polygon to make a cell mask
msk = poly2mask(poly(:,1),poly(:,2),imSize(1),imSize(2));
%index based on points interior to the cell
[r,c] = find(msk == 1);
%measure cell coordinates of each pixel in the cell
[fractionalLength, halfLength, fractionalWidth, halfWidth, ~] = Cell_Projection(cellMesh, [c,r]);

X = fractionalLength*halfLength;
Y = fractionalWidth*halfWidth;

% extract indices of pixels within dC of the cell center line and greater
% than poleLen from the pole.
ix = logical((abs(Y) <= dC).*(X >= min(X)+poleLen).*(X <= max(X)-poleLen));

imageInds = sub2ind(imSize,r(ix),c(ix));
mask = zeros(imSize);
mask(imageInds) = 1;

end