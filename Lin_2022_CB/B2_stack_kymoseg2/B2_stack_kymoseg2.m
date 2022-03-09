
% 2022 Feb 22th: Final version. This code is doing: 
% (1) Use the manually-curated data to refine the linking of separators
% (2) Export kymograph/separator plots for the curated data
% (3) Export a list of separators for manaul curation

% ====================================================================
% Specified the subfolders
SubFolderPrefix = 'waf7_001xy3';
SubFolderList = [1 10]'; % [1 2 3 4 5 6 7 8 9]'; 

% Input parameters. Recommand to use the same parameter as in B1_stack_kymoseg1
param = {};

% kmat construction
param{1}.Nrange = [1 2000];   % Range of frames in image stacks
param{1}.Yrange = [1 425];    % Range of y-axis (rows) in image 

% Filter for good separator.
param{4}.max_frame = param{1}.Nrange(2);     % The last frame number in the data
param{4}.exit_y = param{1}.Yrange(2) - 20;   % threshold of y on exit border
param{4}.minimal_length = 5;                 % minimal length of a separator
param{4}.plot_flag = 1;                      % Use flag==1 for plot
param{4}.plot_separatorA = [1 1];            % Use [1 1] for plot the kymograph/separator (per 1000 frame)
param{4}.plot_separatorB = [1 1];            % Use [1 1] for plot the kymograph/separator (per 100 frame)

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
    [ret1] = sub_kymoseg_step2(path, param);
    
end

% ===================================================================== %
%   Subscript 
% ===================================================================== %

function [ret] = sub_kymoseg_step2(path, param)

cd ( path.save );

% Load previous dataset
temp = load('kymoseg_rec1.mat');
peaks = temp.kymoseg_rec1.peaks;
kmat = temp.kymoseg_rec1.kmat;

% (2-1) Use the manaually-curated .xls to refine the separator linking
[peaks, pair_ind_arrayC, separatorC, peak_linking_tableC] = correct_separator(peaks, path.saveS1, param{4});

% export the plots
[ret2] = export_corrected_separator_data(kmat, peaks, separatorC, param{1}, path.saveS2, param{4}.plot_separatorB);

% (2-2) Export the separator summary for manaul check

[sep_summary] = export_separator_manual_check(separatorC, path.saveS2);

kymoseg_rec2 = {};
kymoseg_rec2.peaks = peaks;
kymoseg_rec2.pair_ind_arrayC = pair_ind_arrayC;
kymoseg_rec2.separatorC = separatorC;
kymoseg_rec2.peak_linking_tableC = peak_linking_tableC;
kymoseg_rec2.sep_summary = sep_summary;
kymoseg_rec2.param = param;
kymoseg_rec2.path = path;

cd(path.save);   save('kymoseg_rec2', 'kymoseg_rec2', '-v7.3');

ret = 1;

end

% ===================================================================== %






