
% Remove debris of segmented region

function [segmatF] = remove_debris(segmat)

segmatF = 0*segmat;


for j = 1: size(segmat, 3)
    
segdata = segmat(:,:,j);

segdataF = 0*segdata;

if (max(max(segdata)) > 0) % having some objects
    
    for m = 1:max(max(segdata))
            
        obj_mask = (segdata == m); 
        
        % Test if this mask has multiple disconnected region
        [mask_output] = keep_largest_object(obj_mask);
                        
        segdataF = segdataF + m * mask_output;
        
    end
    
end

segmatF(:,:,j) = segdataF;

end

%figure;
%subplot(121);  imagesc(segdata);
%subplot(122);  imagesc(segdataF);

end

% =======================================================================

function [mask_output] = keep_largest_object(obj_mask)
      
test_seg = bwlabel(obj_mask);

if (max(max(test_seg)) > 1)  % Having multiple region. Keep the largests
           
    % find the index correspond to the largest region
    
    temp_size = 1;
    index_flag = 1;
    
    for k = 1:max(max(test_seg))
        
        region_size = sum(sum( test_seg == k ));
        
        %region_size
        
        if (region_size > temp_size)
            
            temp_size = region_size;
            index_flag = k;

        end
        
    end
    
    % Only keep the index-k region
    
    mask_output = (test_seg == index_flag);    
    
else
    
    mask_output = obj_mask;
     
end


%figure;
%subplot(121); imagesc(test_seg);
%subplot(122); imagesc(mask_output);

end
