
% 2022 Mar 10th: Final version. This code is doing: 
% Construct "mother cell lineage" (MCL) and append data on the lineage

path = {};
path.folder = 'O:\Shares\Data_04\Wei-Hsiang Lin\uF_dataset\exponential\combined\';
path.save = 'O:\Shares\Data_04\Wei-Hsiang Lin\uF_dataset\exponential\MCL\';
file_name_list = {'M9GlcCA_lac', 'M9GlcCA_ackA', 'M9GlcCA_pta', 'M9GlcCA_adhE', 'M9GlcCA_ldhA'};

param = {};
param.Rsc_adjust_factor = 2;               % Adjust for R_ATP ratio for exposure time (405nm = 200ms, 488 = 400ms)
param.sensor_adjust_factor = [6.84 3.33];  % Weighting factors for QUEEN[2m] protein amount 
param.exposure_adjust_factor = 200;        % exposure time (msec)


% ============================= Script start ===============================

for DB = 1:length(file_name_list)

% (1) Construct MCL (mother cell lineage)

cd(path.folder);
file_name = strcat('combined_', file_name_list{DB});
raw_data = load( strcat(file_name, '.mat')).combined_dataset.cc_ensemble;

lineage_data = {};
flag = 1;

for m = 1:length(raw_data)

    cc_ensemble = raw_data{m};
    
    [CC_data, MCstat, mother_cell] = get_MCL(cc_ensemble);
    
    %%%  Test mother cell lineage %%%
    test1 = ( length(mother_cell.size) > 20 );
    
    if (test1 == 1)
        test2 = ( size(MCstat.cc_data, 1) > 3 );
    end
    
    if ( test1 && test2 )   
        
        [mother_cell] = append_Rsc_sensor(mother_cell, param);
    
        lineage_data{flag}.MCstat = MCstat;
        lineage_data{flag}.MC = mother_cell;    
        
        lineage_data{flag}.param = param;
        lineage_data{flag}.dataset_index = cc_ensemble.dataset_index;
        lineage_data{flag}.dataset_indnum = cc_ensemble.dataset_indnum;
        lineage_data{flag}.chamber_index = cc_ensemble.chamber_index;
        
        flag = flag + 1;   
    end
    
end

% (2) Save the MCL data
cd(path.save);
save_name = strcat('MCL_', file_name);
save( save_name, 'lineage_data' );

end

