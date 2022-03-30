
% 2022 Mar 10th: Final version. This code is doing: 
% Find mother-daughters "triad" structure, which is essential for cell pairs analysis

path = {};
path.folder1 = 'O:\Shares\Data_04\Wei-Hsiang Lin\uF_dataset\exponential\combined\';
path.folder2 = 'O:\Shares\Data_04\Wei-Hsiang Lin\uF_dataset\exponential\CC\';
path.save = 'O:\Shares\Data_04\Wei-Hsiang Lin\uF_dataset\exponential\Tree\';
file_name_list = {'M9GlcCA_lac', 'M9GlcCA_ackA', 'M9GlcCA_pta', 'M9GlcCA_adhE', 'M9GlcCA_ldhA'};

param = {};
param.Rsc_adjust_factor = 2;                % Adjust for R_ATP ratio for exposure time (405nm = 200ms, 488 = 400ms)
param.sensor_adjust_factor = [6.84 3.33];   % Weighting factors for QUEEN[2m] protein amount 
param.exposure_adjust_factor = 200;         % exposure time (msec)

% ============================= Script start ===============================

for DB = 1:length(file_name_list)

% (1) Construct MCL, find tree structure on lineage

cd(path.folder1);
file_name1 = strcat('combined_', file_name_list{DB});
raw_data = load( strcat(file_name1, '.mat')).combined_dataset.cc_ensemble;

cd(path.folder2);
file_name2 = strcat('CC_', file_name_list{DB});
cc_array = load( strcat(file_name2, '.mat')).cc_array;


tree_summary = {};
triad_dataset = {};

for m = 1:length(raw_data)

    cc_ensemble_m = raw_data{m};
    
    [CC_data, MCstat, mother_cell] = get_MCL(cc_ensemble_m);
    
    %%%  Test mother cell lineage %%%
    test1 = ( length(mother_cell.size) > 20 );
    
    if (test1 == 1)
        test2 = ( size(MCstat.cc_data, 1) > 3 );
    end

    if ( test1 && test2 )           
        tree_summary{m} = get_treeV2(CC_data, cc_ensemble_m);          
    else
        tree_summary{m} = NaN;  % tree is too short
    end

end

% (2) From cc_array, append cell division triad data
triad_data = append_triad_data(tree_summary, cc_array);


% (3) Save the tree and triad data

tree = {};
tree.tree_summary = tree_summary;
tree.triad_data = triad_data;

cd(path.save);
save_name = strcat('Tree_', file_name_list{DB});
save( save_name, 'tree' );

end

