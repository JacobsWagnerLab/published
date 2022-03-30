
function [output] = get_treeV2(CC_data, cc_ensemble)

% (1) Construct tree from root
cc_num = size(CC_data,1);

tree_data = NaN(cc_num,3);
% col [1] division degree
% col [2] tree degree
% col [3] aging group

% ----------------------------------------------------------
% (a1) For the first frame, assign division degree as zero
for j = 1:cc_num
    
    if (CC_data(j,1) == 1)  % start from first frame        
        tree_data(j,1) = 0;   % degree zero        
    end
    
end

% ----------------------------------------------------------
% (a2) Extend the division degree
for j = 1:cc_num

    % Find mother cell.    
    mother_id = CC_data(j,3);
    
    % Increase div.degree by 1    
    if (mother_id > 0)        
        tree_data(j,1) =  tree_data(mother_id,1) + 1;
    end

end

% ----------------------------------------------------------
% (a3) Decide the end of MC lineage and trace back
max_div = max(tree_data(:,1));

candidate_list = [];

for j = 1:cc_num
    
    if (tree_data(j,1) == max_div)        
        candidate_list = [candidate_list; j];
    end
    
end

ter_num = length(candidate_list);
yrec = NaN(ter_num,2);  % min, max

for k = 1:ter_num
    
    id = candidate_list(k);
    ydata = cc_ensemble.data{id}.data(:,8);
    yrec(k,:) = [min(ydata) max(ydata)];
    
end

% Find the cell with minimal y
[val,ind] = min(yrec(:,1));  
ter_id = candidate_list(ind);

% Start from terminua MC, trace back
tree_data(ter_id,2) = 1;

temp_id = ter_id;

for k = 1:max_div
    
    mother_id = CC_data(temp_id, 3);
    tree_data(mother_id,2) = 1;
    temp_id = mother_id;
    
end

% ----------------------------------------------------------
%(a4) From MC lineage, append the tree degree 2,3,...etc

max_tree_degree = max_div;

for dg = 1:max_tree_degree
    
   for j = 1:cc_num
    
       if (tree_data(j,2) == dg)
            
           % If there is an offspring not being assgined yet, assign tree degree 
            
           offspring_list = CC_data(j,[4 5]);
           
           for n = 1:2
               
               temp_id = offspring_list(n);
               
               if ( ~isnan(temp_id) )
                   
                   if( isnan(tree_data(temp_id,2)) ) % no assign yet
                        tree_data(temp_id, 2) = dg+1;                               
                   end
                   
               end
               
           end
                      
       end
       
   end
   
end

% ----------------------------------------------------------
%(a5) Generate aging group, write into tree_data 3rd column

for j = 1:cc_num
    
    if (tree_data(j,2) == 1)  % must be age group class 1
        tree_data(j,3) = 1;
    end
    
    if (tree_data(j,2) == 2)  % must be age group class 2
        tree_data(j,3) = 2;
    end
    
end

% Classify age group class (3,4)/(5,6)/(7,8)

for mc_age_class = 2:4  % mother cell age class
    
list = [];

for j = 1:cc_num
    
    % Find mother cell with age group 2. 3. 4
    test1 = (tree_data(j,3) == mc_age_class);    
    test2 = sum( ~isnan(CC_data(j,4:5)) );  % number of daughtors
    
    if (test1 && (test2>0))
        list = [list; j test2 CC_data(j,4:5)];        
    end    
end

list_num = size(list,1);

for k = 1:list_num

    if (list(k,2) == 1)  % Only one daughtor, should be class 3
        if ( ~isnan(list(k,3)) )
            tree_data(list(k,3),3) = (2*mc_age_class)-1;
        elseif ( ~isnan(list(k,4)) )
            tree_data(list(k,4),3) = (2*mc_age_class)-1;
        end
    end
    
    if (list(k,2) == 2)  % Both daughtors exist. Compare y-direction
        
        yA = cc_ensemble.data{list(k,3)}.data(:,8);
        yB = cc_ensemble.data{list(k,4)}.data(:,8);

        yAmin = min(yA);
        yBmin = min(yB);
        
        if (yAmin < yBmin) 
            tree_data(list(k,3),3) = (2*mc_age_class)-1;
            tree_data(list(k,4),3) = 2*mc_age_class;
        elseif (yAmin > yBmin)
            tree_data(list(k,3),3) = 2*mc_age_class;
            tree_data(list(k,4),3) = (2*mc_age_class)-1;          
        else
            tree_data(list(k,3),3) = NaN;
            tree_data(list(k,4),3) = NaN; 
        end
        
    end
    
end

end

% ----------------------------------------------------------
% (b) Generate triad list for mother-offspring data
% Note that this list is only for this chamber.

triad_list = [];
% [1-3] mother, offspring1, offspring2, 
% [4-6] age class of each cells

for j = 1:cc_num
    
    offspring_list = CC_data(j,[4 5]);    
    test = sum(~isnan(offspring_list));

    if ( test == 2 )  % Both are not NaN
        
        write = [j offspring_list];
        triad_list = [triad_list; write];
    end

end

for r = 1:size(triad_list)
    
    triad_list(r,4) = tree_data(triad_list(r,1), 3);
    triad_list(r,5) = tree_data(triad_list(r,2), 3);
    triad_list(r,6) = tree_data(triad_list(r,3), 3);
    
end

%==========================================================

output = {};
output.CC_data = CC_data;
output.tree_data = tree_data;
output.triad_list = triad_list;

output.dataset_index = cc_ensemble.dataset_index;
output.dataset_indnum = cc_ensemble.dataset_indnum;
output.chamber_index = cc_ensemble.chamber_index;

end



