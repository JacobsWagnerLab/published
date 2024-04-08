
% Find level_matrix of an object
% input_mat may has binray data

function [level_mat] = find_level_mat(input_mat, max_level)

% Convert into double matrix
input_mat = double(input_mat);
level_mat = input_mat;

for r = 1:size(input_mat,1)
    for c = 1:size(input_mat,2)
        
        if (input_mat(r,c) > 0)
            level_mat(r,c) = NaN;
        end        
    end
end

for g = 1:max_level
       
    for r = 1:size(level_mat,1)
        
        for c = 1:size(level_mat,2)
            
            if ( isnan(level_mat(r,c)) == 1)  % unassigned, looking for its neighbor
                
                if (r == 1)                       % upper boundary
                   
                    level_mat(r,c) = g;
                
                elseif (r == size(level_mat,1))   % lower boundary
                    
                    level_mat(r,c) = g;          
                    
                elseif (c == 1)                   % left boundary
                    
                    level_mat(r,c) = g;         
                    
                elseif (c == size(level_mat,2))   % right boundary
                    
                    level_mat(r,c) = g;
                
                else  % Not a boundary pixel, search for neighborhood
                                        
                    if (level_mat(r-1,c) == (g-1) )
                        
                        level_mat(r,c) = g;
                    
                    elseif (level_mat(r+1,c) == (g-1) )
                        
                         level_mat(r,c) = g;
                         
                    elseif (level_mat(r,c-1) == (g-1) )    
                        
                         level_mat(r,c) = g;
                    
                    elseif (level_mat(r,c+1) == (g-1) ) 
                         
                         level_mat(r,c) = g;

                    end
                        
                end
                
            end
            
        end
        
    end
    
end


