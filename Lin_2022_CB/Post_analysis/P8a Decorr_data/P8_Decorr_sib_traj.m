
% Assuembly the ATP trajecotry data between siblings
% Exported as "ATP_traj_data".

path.load = '/Volumes/Data_04/Wei-Hsiang Lin/uF_dataset/exponential/Tree/';

dataset_name = {'M9GlcCA_lac', 'M9GlcCA_ldhA', 'M9GlcCA_adhE', 'M9GlcCA_pta', 'M9GlcCA_ackA', ...
                'M9Glc_lac', 'M9Xyl_lac', 'M9Gly_lac'};

% ======================= Script starts ================================

ATP_traj_data = {};

for DB = 1:length(dataset_name)

    cd(path.load);
    
    tree_data = load( strcat('Tree_', dataset_name{DB}) ).tree;
    triad0 = tree_data.triad_data;
    
    ATP_traj_data{DB}.data = get_ATP_traj_siblings(triad0);    
    ATP_traj_data{DB}.dataset_name = dataset_name{DB};
        
end


% =====================================================================

function [output] = get_ATP_traj_siblings(triad0)

% (1) Filter for parent and both siblings to be present
[triad0] = filter_triad(triad0, 9, 1);
[triad0] = filter_triad(triad0, 10, 1);
[triad0] = filter_triad(triad0, 11, 1);

% (2) Extract triade data of ATP trajectory

ATP_traj = {};
type_flag = ones(4,1);

triad_num = size(triad0.rec,1);

for k = 1:triad_num

    mtype = triad0.rec(k,4);
    
    if ( triad0.rec(k,6) > triad0.rec(k,5) )  % ccB, ccC
        
        ATP_traj{mtype}{type_flag(mtype)}.ccP = triad0.triad_dataset{k}.ccA.Rsc_traj;
        ATP_traj{mtype}{type_flag(mtype)}.cc1 = triad0.triad_dataset{k}.ccB.Rsc_traj;
        ATP_traj{mtype}{type_flag(mtype)}.cc2 = triad0.triad_dataset{k}.ccC.Rsc_traj;        
        
        type_flag(mtype) = type_flag(mtype)+1;
        
    elseif (  triad0.rec(k,6) < triad0.rec(k,5)  )  % ccC, ccB
        
        ATP_traj{mtype}{type_flag(mtype)}.ccP = triad0.triad_dataset{k}.ccA.Rsc_traj;
        ATP_traj{mtype}{type_flag(mtype)}.cc1 = triad0.triad_dataset{k}.ccC.Rsc_traj;
        ATP_traj{mtype}{type_flag(mtype)}.cc2 = triad0.triad_dataset{k}.ccB.Rsc_traj; 
        
        type_flag(mtype) = type_flag(mtype)+1;
        
    end
    
end

% (3) Calcualte difference of ATP trajs between sibling pairs

max_type = length(ATP_traj);

for mtype = 1:max_type

    for j = 1:length(ATP_traj{mtype})
        
        data1 = ATP_traj{mtype}{j}.cc1;
        data2 = ATP_traj{mtype}{j}.cc2;
        L1 = length(data1);
        L2 = length(data2);
        Lmin = min(L1,L2);
        
        ATP_traj{mtype}{j}.diff = data2(1:Lmin) - data1(1:Lmin);        
        
    end
    
end

% (4) Collect ATP trajectory difference between siblings into matrix 

ATP_traj_mat = {};
ATP_traj_mat2 = {}; % absolute difference 

max_length = 200;

for mtype = 1:max_type
    
    ATP_traj_mat{mtype} = NaN( length(ATP_traj{mtype}), max_length);
    ATP_traj_mat2{mtype} = NaN( length(ATP_traj{mtype}), max_length);
    
    for j = 1:length(ATP_traj{mtype})
    
        L = length(ATP_traj{mtype}{j}.diff);
        ATP_traj_mat{mtype}(j,1:L) = ATP_traj{mtype}{j}.diff;
        ATP_traj_mat2{mtype}(j,1:L) = abs(ATP_traj_mat{mtype}(j,1:L));
        
    end
    
end

% (5) Collect ATP trajectory difference 
%    (Combine all different types, and polarized by {older cell}- {newer cell})

ATP_traj_mat_total = [];
ATP_traj_mat2_total = [];

for mtype = 1:2

    if (mtype == 1)
        
        write = ATP_traj_mat{mtype}; 
        write2 = ATP_traj_mat2{mtype};
        %size(write)
        %size(write2)
        ATP_traj_mat_total = [ATP_traj_mat_total; write];
        ATP_traj_mat2_total = [ATP_traj_mat2_total; write2];
        
    elseif (mtype == 2)
        
        write = (-1)* ATP_traj_mat{mtype}; 
        write2 = ATP_traj_mat2{mtype};
        ATP_traj_mat_total = [ATP_traj_mat_total; write];
        ATP_traj_mat2_total = [ATP_traj_mat2_total; write2];
        
    end
    
end

% (5) Save data

output = {};
output.ATP_traj = ATP_traj;
output.ATP_traj_mat = ATP_traj_mat;
output.ATP_traj_mat2 = ATP_traj_mat2;
output.ATP_traj_mat_total = ATP_traj_mat_total;
output.ATP_traj_mat2_total = ATP_traj_mat2_total;


end

% =====================================================================

function [triadN] = filter_triad(triad, filter_col, filter_label)

flag = 1;
triad_num = size(triad.rec,1);

triadN = {};
triadN.triad_dataset = {};
triadN.rec = [];

for j = 1:triad_num
    
    rtemp = triad.rec(j,:);
    
    if (rtemp(filter_col) == filter_label) % find matched item
        
        triadN.triad_dataset{flag} = triad.triad_dataset{j};
        triadN.rec(flag,:) = triad.rec(j,:);        
        flag = flag + 1;
        
    end
    
end

end

% =====================================================================
