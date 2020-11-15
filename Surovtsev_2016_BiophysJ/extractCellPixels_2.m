function [roi_indx, px_values,roi] = extractCellPixels_2(cellImage,cellMesh, dims)
%--------------------------------------------------------------------------
% returns linear indices of the cellImage which falls into region of
% interest, ROI
% ROI is defined relative to centerline in the cellMesh - at given location 
% along the centerline and within given radius 
%--------------------------------------------------------------------------
%
%**********INPUT********:
% cellImage = image of the cell (i.e. one frame aquired at the microscope)
% cellMesh = cell mesh coordinates as defined in oufti's cellList (should correspond to the cellImage) 
% dims = 3-value array specifiying ROI:
%       dims(1) and dims(2) are initial and final coordinates along the centerline;
%       dims(3) is a radius of the ROI from the centerline 
%
%*********OUTPUT********:
% roi_indx: linear indices of the pixels within ROI in cellImage.
% px_values: values of the pixels within ROI in cellImage.
% roi: coordinates (in image coordinates) of ROI
%
%@author:  Ivan Surovtsev
%@date:    July 29, 2015
%@copyright 2012-2015 Yale University
%

im_size=size(cellImage);
% centerline coordinates
cline=[mean(cellMesh(:,[1,3]),2),mean(cellMesh(:,[2,4]),2)];
 r=dims(3); % radius around centerline

% to calculate norm of the vectors
dist=@(xy1,xy2) sqrt(sum((xy2-xy1).^2)); 

% lets get ROI
% *01, *02 - refers to start end points on centerline
% *11,*12 - refers to edge points of the ROI
% *m* - refers to corresponidng vetices in meshes
xy_01=cline(dims(1),:); 
 xy_m11=cellMesh(dims(1),[1,2]);
  xy_m12=cellMesh(dims(1),[3,4]);
xy_02=cline(dims(2),:);
 xy_m21=cellMesh(dims(2),[1,2]);
  xy_m22=cellMesh(dims(2),[3,4]);

xy_11=xy_01-r*(xy_m12-xy_m11)/dist(xy_m12,xy_m11);
 xy_12=xy_01+r*(xy_m12-xy_m11)/dist(xy_m12,xy_m11);
xy_21=xy_02-r*(xy_m22-xy_m21)/dist(xy_m22,xy_m21);
 xy_22=xy_02+r*(xy_m22-xy_m21)/dist(xy_m22,xy_m21);

roi=double([xy_11;xy_12;xy_22;xy_21]);

% make a mask over ROI
msk = poly2mask(roi(:,1),roi(:,2),im_size(1),im_size(2));
 
% getting indices and values
lin_indx=1:length(cellImage(:));
 roi_indx=lin_indx(msk==1);
px_values=double(cellImage(msk==1));

end