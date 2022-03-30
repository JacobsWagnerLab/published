
function [mask] = OJ_reindex(input)

% Re-index object such that the each object index is nonempty
% input: a labeled matrix
% ouput: a labeled matrix with filtered objects

obj_num = max(max(input));
record = zeros(obj_num, 2);

% find the object size for each index

for m = 1:obj_num    
    record(m) = sum(sum(input == m));        
end

flag = 0;

% setup new index

for m = 1:obj_num
    
    if ( record(m) > 0 )
        flag = flag + 1;
        record(m,2) = flag;        
    end
    
end

% draw the labeled matrix with new index

mask = 0*input;

for m = 1:obj_num
    
    mask_m = (input == m);     
    new_index = record(m,2);
    
    if (new_index > 0)
    
        mask = mask + mask_m * new_index;
    
    end
    
end

%  
%  figure;
%  subplot(121); imagesc(input);
%  subplot(122); imagesc(mask);
%  colormap(jet); 
% 


