
function [output] = RGB_format(input, colormap)

% Generate RGB matrix for figure matlab imwrite() function
% colormap = hot;
% input = C1;

%%%

cM = colormap; 
img_size = size(input);
output = zeros(img_size(1), img_size(2), 3);

input = double(input);
A_norm = input / ceil(max(max(input)));         
    
for r = 1: img_size(1)
        
    for c = 1: img_size(2)
            
        for color = 1:3
                
            ind_temp = floor(256 * A_norm(r,c)); 
            ind_temp = max(ind_temp,1);    
            output(r,c,color) = cM(ind_temp, color);  
            
        end   
    end
end
    
%imagesc(output)




