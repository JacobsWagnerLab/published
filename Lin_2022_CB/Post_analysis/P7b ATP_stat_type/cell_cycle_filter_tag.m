
% Filter particular tag for cell cycle array;

function [cc_arrayN] = cell_cycle_filter_tag(cc_array, tagA, tagB, range)

% Example
%tagA = 'GR_stat';
%tagB = 'GR_CoD';
%range = [0.7 1];

cc_arrayN = {};
flag = 1;

for j = 1:length(cc_array)
    
    if ( ~isempty(cc_array{j}.(tagA)) )
    
        temp = cc_array{j}.(tagA).(tagB);
    
        if ( ~isempty(temp) )
    
            if ( temp > range(1) ) && ( temp < range(2) )
        
                cc_arrayN{flag} = cc_array{j};
                flag = flag + 1;
                
            end
            
        end
        
    end
    
end
