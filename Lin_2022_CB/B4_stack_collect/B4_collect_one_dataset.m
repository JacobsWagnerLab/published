
% 2022 Feb 22th: Final version. This code is doing: 
% Combine all cell cycle data of all chambers in single field of view (FOV).

% Input data path information:

path.server = 'O:/Shares/Data_04/Wei-Hsiang Lin/WHLin_data_N/';
path.folder = '/KymoSeg/';
path.save = 'O:/Shares/Data_04/Wei-Hsiang Lin/uF_dataset/exponential/FOV/';

save_nametag = '20220114B';
folder_name = '20220114/waf7_002xy2-';
subfolder = [1:5];

% ====================================================================== %
% Scrip start
% ====================================================================== %

cc_ensemble = {};

chamber_flag = 1;

for j = 1:length(subfolder)
    
    target_folder = strcat(path.server, folder_name, num2str(subfolder(j)), path.folder);

    cd(target_folder);
    temp = load('kymoseg_rec3.mat').kymoseg_rec3;   
    max_frame = max( temp.obj_link_list(:,1) );
    
    cc_ensemble{chamber_flag}.data = temp.cc_ensemble;
    cc_ensemble{chamber_flag}.max_frame = max_frame;
    cc_ensemble{chamber_flag}.track_note = temp.cc_note;
    cc_ensemble{chamber_flag}.track_check = temp.cc_check;
    cc_ensemble{chamber_flag}.chamber_index = subfolder(j);
    cc_ensemble{chamber_flag}.dataset_index = save_nametag;
    
    chamber_flag = chamber_flag + 1;
    
end

cd(path.save);
save_name = strcat('cc_ensemble_', save_nametag);
save(save_name, 'cc_ensemble', '-v7.3');

%=====================================================================

