
function [mask_ind] = Mask_terra(input, cap_param)

% Making territary mask, need input a segmented matrix

%input = B3(100:200, 100:200);
%input=  B3;
%cap_param = [7 5];

seg_core = (input>0);
core_label = bwlabel(seg_core);
cap_terra = PM_2D_gaussian(cap_param(1), cap_param(2));
obj_num = max(max(core_label));

%

mask_ind = 0 *core_label;
mask_val = 0 *core_label;

for j = 1: obj_num
    
    coreJ = double(core_label == j);    
    terraJ = ( conv2(coreJ, cap_terra, 'same') );
    
    for r = 1:size(seg_core,1)
    
        for c = 1:size(seg_core,2)
            
            val_temp = terraJ(r,c);
            
            if ( val_temp > mask_val(r,c) )   
                
                % find a greater candidate for terretary
                % replace the affilation of this pixel to object j
                
                mask_val(r,c) = val_temp;
                mask_ind(r,c) = j;
                                
            end
            
            
        end
        
    end
        
end

%%% 


plotflag = 0;

if (plotflag == 1)

    figure;

    subplot(231);   imagesc(input);
    subplot(232);   imagesc(core_label);
    subplot(233);   imagesc(mask_ind);
    subplot(234);   imagesc(mask_val);

    colormap(jet);

end


