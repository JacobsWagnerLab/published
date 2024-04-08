
function [mask] = OJ_size_filter(input, size_filter)

% Filter small or large object
% input: a labeled matrix
% ouput: a labeled matrix with filtered objects

%input = C2;
%size_filter = [80 1000];  % lower bound, upper bound

mask = input;
obj_num = max(max(mask));

for m = 1:obj_num
    
    obj_m = (mask == m);   % binary mask for object_m
    obj_size = sum(sum(obj_m));

    flag = ( (obj_size < size_filter(1)) || (obj_size > size_filter(2)) );
    
    if (flag == 1) % this object is either too small or too large
        
        mask = mask .* (obj_m == 0);
        
    end    
    
end

 
% figure;
% subplot(121); imagesc(input);
% subplot(122); imagesc(mask);
% 

