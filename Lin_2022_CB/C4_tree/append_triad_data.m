
function [output] = append_triad_data(tree_summary, cc_array)

% Append triad dataset on triad_dataset

CB_num = length(tree_summary);

triad_dataset = {};
flag = 0;
rec = [];
% [1-3] mother, daughtorA, daughtorB index (in each chamber)
% [4]   division degree
% [5]   tree degree
% [6]   aging group
% [7]   FOV number
% [8]   chamber number index (note: not CB)
% [9-11]  whether the cc_array data is found for mother, daughtorA, daughtorB

for m = 1:CB_num
    
if ( isfield(tree_summary{m}, 'CC_data' ) )
    
dataset_indnum = tree_summary{m}.dataset_indnum;
chamber_index = tree_summary{m}.chamber_index;
triad_list = tree_summary{m}.triad_list;

cc_num = size(triad_list,1);
serach_result = zeros(cc_num,3);

% ======================================================================

for j = 1:cc_num   
        
    flag = flag + 1;
    triad_dataset{flag}.ccA = NaN;
    triad_dataset{flag}.ccB = NaN;
    triad_dataset{flag}.ccC = NaN;  
            
    % Search for the cc_array, see if the cell cycle data is present
    
    searchA = [dataset_indnum chamber_index triad_list(j,1)];
    searchB = [dataset_indnum chamber_index triad_list(j,2)];
    searchC = [dataset_indnum chamber_index triad_list(j,3)];
    
    for cc = 1:length(cc_array)
            
        temp = cc_array{cc}.cc_info;
        target = [temp.FOV temp.chamber_index cc_array{cc}.ID(3)];
        
        if ( searchA == target )
            triad_dataset{flag}.ccA = cc_array{cc};
            serach_result(j,1) = 1;
        end
        
        if ( searchB == target )
            triad_dataset{flag}.ccB = cc_array{cc};
            serach_result(j,2) = 1;
        end        
        
        if ( searchC == target )
            triad_dataset{flag}.ccC = cc_array{cc};
            serach_result(j,3) = 1;
        end
        
    end
    

end

tempvec = ones(cc_num,1);
write = [triad_list tempvec*dataset_indnum tempvec*chamber_index serach_result];
rec = [rec; write];

% ======================================================================

end

end

output.triad_dataset = triad_dataset;
output.rec = rec;

end

