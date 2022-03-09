
function [mask_output] = Mask_fill(mask)

% Fill the inner holes of the object on the mask
%mask = mask_record.SegMatContent(:,:,10);

% (1) Invert 0/1 of the mask
maskR = (mask == 0);

% (2) Finding the "background component" of the maskR
segmatR = bwlabel(maskR);
bg_index = segmatR(1,1);  % Using the upper left corner 

mask_output = mask;
mask_output( (segmatR ~= bg_index) ) = 1;

%figure;  
%subplot(121); imagesc(mask);
%subplot(122); imagesc(mask_output);