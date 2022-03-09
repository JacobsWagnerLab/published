
% 2022 Feb 22th: Final version. This code is doing: 
% (1) Used the refined separator data to aid cell segmentation
% (2) Export the cell cycle data and append fluorescence data

% ====================================================================
% Specified the subfolders
SubFolderPrefix = 'waf7_001xy3';
SubFolderList = [1 10]'; % [1 2 3 4 5 6 7 8 9]'; 

% Input parameters
param = {};
param{5}.sep_delay = 20;        % Delay frame of separator
param{5}.frame_interval = 6;    % Frame interval
param{5}.max_digit = 4;         % Maximal digit in dataset
param{5}.append_fluo = 1;       % Use append_fluo==1 for adding fluorescence data

name = {};      
name.c2.prefix = strcat(SubFolderPrefix ,'t');
name.c2.postfix = 'c2s.tif';
name.c3.prefix = strcat(SubFolderPrefix ,'t');
name.c3.postfix = 'c3s.tif';
param{5}.name = name; 

% ====================================================================

close all;

parfor j = 1:size(SubFolderList,1)
       
    TagJ = num2str( SubFolderList(j));
    PathJ = strcat( data_folder, SubFolderPrefix, '-', TagJ, '\' );
    
    path = {};
    path.load = strcat(PathJ, 'PreSegB\');
    path.c2 = strcat(PathJ, 'c2Exp\');
    path.c3 = strcat(PathJ, 'c3Exp\');    
    path.save = strcat(PathJ, 'KymoSeg\');
    path.saveS1 = strcat(PathJ, 'KymoSeg\Separator\');
    path.saveS2 = strcat(PathJ, 'KymoSeg\SeparatorC\');
    path.saveSeg = strcat(PathJ, 'KymoSeg\SegN\');
    
    cd(script_dir);       
    [ret1] = sub_kymoseg_step3(path, param);
    
end


% ===================================================================== %
%   Subscript 
% ===================================================================== %

function [ret] = sub_kymoseg_step3(path, param)

cd(path.load);
temp = load('mask_record.mat');
segmatD = temp.mask_record.SegMatContentD;   % cell outline

cd(path.save);  
temp = load('kymoseg_rec2.mat');
peaks = temp.kymoseg_rec2.peaks;    
separatorC = temp.kymoseg_rec2.separatorC;    
sep_summary = temp.kymoseg_rec2.sep_summary; 
param1 = temp.kymoseg_rec2.param{1};
param4 = temp.kymoseg_rec2.param{4};

% (3-1) Convert into kseg matrix
 [kseg, sep_filter] = separator_to_kseg_V2(peaks, sep_summary, param1, param{5}, path.saveS2);

% (3-2) Print and save the new segmented stacks
[segmatN] = kseg_to_segstack(segmatD, kseg, param1);

if (param4.print_segmat == 1)
    [retP] = print_segmat_panel(segmatN, path.saveSeg, param1);
end

% (3-3) Construct cell cycle based on separator
[obj_link_list, cc_note, cc_check, cc_ensemble] = construct_cell_cycle(separatorC, sep_filter, segmatN); 

% (3-4) Append c2, c3 intensity

if (param{5}.append_fluo == 1)
    
    fluo_rec = measure_intensity(segmatN, path, param1, param{5});    
    [cc_ensemble] = append_cc_intensity(cc_ensemble, fluo_rec);

end

[ret1] = plot_cell_cycle(cc_ensemble, cc_check, path.save);

% Save the intermediate result for second step.
kymoseg_rec3 = {};
kymoseg_rec3.kseg = kseg;
kymoseg_rec3.obj_link_list = obj_link_list;
kymoseg_rec3.cc_note = cc_note;
kymoseg_rec3.cc_check = cc_check;
kymoseg_rec3.cc_ensemble = cc_ensemble;
kymoseg_rec3.param = param;
kymoseg_rec3.path = path;

cd(path.save);   save('kymoseg_rec3', 'kymoseg_rec3', '-v7.3');

% Append segmatN on the mask dataset
mask_record.SegMatContentN = segmatN;
cd(path.save);   save('mask_recordN', 'mask_record', '-v7.3');

%}

ret = 1;

end

% ===================================================================== %






