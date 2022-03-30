
% This code collect all "triad" (mother cell and daughter cell pairs)
% in the tree structure and generate data structure "triad_stat"
% for following analysis

path.load = '/Volumes/Data_04/Wei-Hsiang Lin/uF_dataset/exponential/Tree/';

dataset_name = {'M9GlcCA_lac', 'M9GlcCA_ldhA', 'M9GlcCA_adhE', 'M9GlcCA_pta', 'M9GlcCA_ackA', ...
                'M9Glc_lac', 'M9Xyl_lac', 'M9Gly_lac'};

triad_stat = {};

tag_array = {};
tag_array{1} = {'GR_stat', 'GR'};
tag_array{2} = {'Rsc_stat', 'mean'};
tag_array{3} = {'Rsc_stat', 'std'};
tag_array{4} = {'Rsc_stat', 'ini'};
tag_array{5} = {'Rsc_stat', 'end'};
tag_array{6} = {'size_stat', 'ini'};
tag_array{7} = {'size_stat', 'add1'};
tag_array{8} = {'sensor_stat', 'mean'};


% ======================= Script starts ================================

for DB = 1:length(dataset_name)

    cd(path.load);
    
    tree_data = load( strcat('Tree_', dataset_name{DB}) ).tree;
    triad0 = tree_data.triad_data;
    
    triad_stat{DB}.data = triad_statistics(triad0, tag_array);
    triad_stat{DB}.dataset_name = dataset_name{DB};
    
    
end


% =====================================================================

function [triad_array] = triad_statistics(triad0, tag_array)

% (1) Filter for both silbings to be present
[triad0] = filter_triad(triad0, 9, 1);
[triad0] = filter_triad(triad0, 10, 1);
[triad0] = filter_triad(triad0, 11, 1);

% (2) Separate the triad data according to mother age class 1,2,3,4
triad_array = {};

for k = 1:4
    triad_array.data{k} = filter_triad(triad0, 4, k);
end

tag_type = length(tag_array);

for  t = 1:tag_type

    tag1 = tag_array{t}{1};
    tag2 = tag_array{t}{2};
    
    for k = 1:4
        triad_array.(tag1).(tag2){k} =  extract_triad_data(triad_array.data{k}, tag1, tag2);    
    end
    
end

end

% =====================================================================

function [output] = extract_triad_data(data, tag1, tag2)

SBL_num = size(data.rec,1);  % number of siblings pairs

output = NaN(SBL_num, 4);  %[1,2] data;  [3,4] sibling age class

for j = 1:SBL_num
    
    output(j,1) = data.triad_dataset{j}.ccA.(tag1).(tag2);
    output(j,2) = data.triad_dataset{j}.ccB.(tag1).(tag2);
    output(j,3) = data.triad_dataset{j}.ccC.(tag1).(tag2);
    output(j,4) = data.rec(j,4);
    output(j,5) = data.rec(j,5);
    output(j,6) = data.rec(j,6);
    
end

% Sort the data with increaseing age group

for j = 1:SBL_num
    
    if (output(j,5) > output(j,6))  % revert the data order
        
        temp = [output(j,1) output(j,3) output(j,2) output(j,4) output(j,6) output(j,5)];
        output(j,:) = temp;
    end
        
end

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
