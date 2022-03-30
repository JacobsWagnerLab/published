
% Classify cell cycle according to age class (type 1 to 8)
% Manaully load sibling_stat.m

cc_array_type = {};

for DB = 1:8

ccdata = sibling_stat{DB}.data.data;

cc_arrayT = {};          % type 1 to type 8
count_flag = zeros(8,1); % counting flag from type 1 to type 8


for mtype = 1:4  % mother type of 1 to 4 
    
    rec = ccdata{mtype}.rec;
    temp = ccdata{mtype}.triad_dataset;
    num = length(temp);
    
    for k = 1:num  % go through all triads
        
        % Checking the rec to see which the cell type of ccB and ccC
        
        % Get type of (2m-1) cell cycle data from ccB, 
        % and type of (2m) cell cycle data from ccC.
        
        type_ccB = rec(k,5);
        type_ccC = rec(k,6);

        count_flag(type_ccB) = count_flag(type_ccB) + 1;
        cc_arrayT{type_ccB}{count_flag(type_ccB)} = temp{k}.ccB;
        
        count_flag(type_ccC) = count_flag(type_ccC) + 1;
        cc_arrayT{type_ccC}{count_flag(type_ccC)} = temp{k}.ccC;
         
    end
    
    
end

cc_array_type{DB}.cc_array = cc_arrayT;
cc_array_type{DB}.dataset_name = sibling_stat{DB}.dataset_name;

end


