
% impose upper/lower limit of an image
% used for pre-threshold the inclusion body 

function [Img] = Img_LU_cut(Img0, CutRg)

Img = Img0;

for j1 = 1: size(Img,1)
    
    for j2 = 1:size(Img,2)
                
        if ( Img(j1,j2) > CutRg(2) )    % brighter than upper limit
            
            Img(j1,j2) = CutRg(2);
            
        elseif ( Img(j1,j2) < CutRg(1) )    % dimmer than lower limit
          
            Img(j1,j2) = CutRg(1);
        
        end
        
    end
    
end



