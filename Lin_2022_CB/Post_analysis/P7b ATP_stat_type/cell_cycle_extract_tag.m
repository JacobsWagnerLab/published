
% Extract given quantity of cell cycle array

function [output] = cell_cycle_extract_tag(cc_array, tagA, tagB)

% Example:
%tagA = 'Rsc_stat';
%tagB = 'mean';


output = NaN(length(cc_array), 1);

for j = 1:length(cc_array)
    
    if ( ~isempty(cc_array{j}.(tagA)) )
    
        temp = cc_array{j}.(tagA).(tagB);
    
        if ( ~isempty(temp) ) 
                        
            output(j) = temp;
            
        end
        
    end
    
end
