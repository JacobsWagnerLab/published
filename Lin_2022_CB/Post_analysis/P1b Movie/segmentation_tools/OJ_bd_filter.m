
function [mask] = OJ_bd_filter(input, bd_param)

% Filter small or large object
% input: a segmented matrix
% ouput: a labeled, segmented matrix with filtered objects

%input = B3;
%bd_width = 2;   % thickness of boundarys
%OL_size = 3;    % maximal pixel overlapped with boundary region

bd_width = bd_param(1);
OL_size = bd_param(2);

mask0 =  input; 
img_size = size(mask0);

mask = mask0;
obj_num = max(max(mask));

bd_mask = 0 * mask0;
bd_mask(1:bd_width, :) = 1;
bd_mask(img_size(1)-bd_width+1 : img_size(1), :) = 1;
bd_mask(:, 1:bd_width) = 1;
bd_mask(:, img_size(2)-bd_width+1 : img_size(2) ) = 1;

for m = 1:obj_num
    
    obj_mask = double(mask == m);
    
    flag = sum(sum(bd_mask .* obj_mask));
    
    if (flag > OL_size ) % this object is overlapped with boundary

        mask = mask .* (obj_mask == 0);
        
    end  
    
end

% 
% figure;
% subplot(121); imagesc(mask0);
% subplot(122); imagesc(mask);


