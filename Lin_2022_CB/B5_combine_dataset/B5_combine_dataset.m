
% 2022 Feb 22th: Final version. This code is doing: 
% Combine all cell cycle data of all FOV with same experimental condition

% Input data path information:

path.server = '/Shares/Data_04/Wei-Hsiang Lin/WHLin_data_N/';
path.save1 = '/Shares/Data_04/Wei-Hsiang Lin/uF_dataset/exponential/FOV/';
path.save2 = '/Shares/Data_04/Wei-Hsiang Lin/uF_dataset/exponential/combined/';

combined_name = 'M9GlcCA_ackA_FO';

name_tag = {};
name_tag{1} = '20210423A';
name_tag{2} = '20210423B';
name_tag{3} = '20210423C';
name_tag{4} = '20210423D';
name_tag{5} = '20210423E';
name_tag{5} = '20210423F';

% ========================================================================

name_prefix = 'cc_ensemble_';

combined_dataset = {};
cc_ensemble = {};

cd(path.save1); 
 
flag = 1;

for j = 1:length(name_tag)
          
    dataset_name = strcat(name_prefix, name_tag{j});
    dataset_j = load(dataset_name).cc_ensemble;
    
    for k = 1:length(dataset_j)    
        cc_ensemble{flag} = dataset_j{k};    
        cc_ensemble{flag}.dataset_indnum = j;
        flag = flag + 1;  
    end    
    
end

combined_dataset.cc_ensemble = cc_ensemble;
combined_dataset.tag_list = name_tag;
combined_dataset.name = combined_name;

cd(path.save2);

save_name = strcat('combined_', combined_name, '.mat');
save(save_name, 'combined_dataset');

