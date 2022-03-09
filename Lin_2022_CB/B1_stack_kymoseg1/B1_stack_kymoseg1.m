
% 2022 Feb 21th: Final version. This code is doing: 
% (1) Create kymograph for segmentation (see document for the detail)
% (2) Find "Separators" on kymograph 
% (3) Export kymograph/separator plots and Excel files for manual curation

% ====================================================================

% Specified the subfolders
SubFolderPrefix = 'waf7_001xy3';
SubFolderList = []';

% Input parameters
param = {};

% kmat construction
param{1}.Nrange = [1 2000];  % Range of frames in image stacks
param{1}.Yrange = [1 425];   % Range of y-axis (rows) in image 

% 1D convolution for kmat
param{2}.L = 25;             % Radius for Gaussian convolution
param{2}.sig = 3;            % SD     for Gaussian convolution

% Finding peaks
param{3}.Nrange = param{1}.Nrange;   % Column range of kmat; equal to range of frames 
param{3}.Nsubrange = [1 2000];       % Column range for peak-finding; can be a subset 
param{3}.append_option = 1;          % Overwrite option. Set equal to 1 for overwriting a subset of peak data
param{3}.min_level = (-8);           % Minimal level for peak finding       
param{3}.w = 10;                     % Radius (in pixel) for peak finding 
param{3}.offset = 0.7;               % Peak threshold parameter. Peak must larger than (mean + offset*SD )

% Link peaks between frames, with a compound distance metric 
param{3}.dist_cut = 20;              % minimal distance for peak pairs      
param{3}.dist_weight = [1 1];        % weighting factor for [ypos, peak intensity]
param{3}.dist_cutW = 40;             % minimal distance for weak peak pair

% Filter for good separator
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
    
    cd (PathJ);     mkdir('KymoSeg');
    cd (path.save); mkdir('Separator');
    cd (path.save); mkdir('SeparatorC');
    cd (path.save); mkdir('SegN');
    
    cd(script_dir);       
    [ret1] = sub_kymoseg(path, param);
    
end


% ===================================================================== %
%   Subscript 
% ===================================================================== %

function [ret] = sub_kymoseg(path, param)
    
cd (path.load);
temp = load('mask_record.mat');
segmat = temp.mask_record.SegMatContent;  % Segmentation data of "core region"

% (0) Create 2D kymograph segmentation matrix

[kmat0] = stack_to_kymo(segmat, param{1});
kmat0 = (-1)*kmat0;  % invert the 2D kymomatrix


% (1) Convolve the kymo-matrix with 1D gaussian kernel

kmat_ini = kmat0(:,param{1}.Nrange(1):param{1}.Nrange(2));
kmat = 0*kmat_ini;

for c = 1:size(kmat,2)
    
    ctemp = kmat_ini(:,c);
    gK = PM_1D_gaussian(param{2}.L, param{2}.sig);    
    kmat(:,c) = conv(ctemp, gK, 'same');

end

% (2) Perform peak detection on kmat. 
%     This code supports partial-overwrite modification.
%     User can set new parameter for a subset of frames and find peaks again.

cd(strcat(path.save, '\Separator')); 

if ( ~isfile('peaks_data.mat') )
    disp('create new variable');
    peak = {};
    for fn = param{1}.Nrange(1):param{1}.Nrange(2)
        peaks{fn}.data = [];
    end
    save('peaks_data.mat', 'peaks');
end

peaks = load('peaks_data.mat').peaks;

[peaks] = kmat_to_peaks(kmat, peaks, param{3});   
[pair_ind_array, peaks] = connect_peaks(peaks, param{3});

% (3) Link peaks into separators
[peaks, separator] = peak_to_separator(peaks, path.saveS1, param{4});
[peak_linking_table] = export_separator_data(kmat, peaks, param{1}, separator, path.saveS1, param{4}.plot_separatorA);


% (4) Save the intermediate result for the first step

% (4-1) Save the parameter for each frame
param_array = {};

if ( isfile('kymoseg_rec1.mat') )  % Loading param array data (if existed)  
    param_array = load('kymoseg_rec1.mat').kymoseg_rec1.param_array;
else
    % Do nothing
end

for j = param{3}.Nsubrange(1):param{3}.Nsubrange(2)
    param_array{j} = param{3};    
end

cd(path.saveS1); 
save('peaks_data.mat', 'peaks');

% (4-2) Save all data

kymoseg_rec1 = {};
kymoseg_rec1.segmat = segmat;
kymoseg_rec1.kmat = kmat;
kymoseg_rec1.peaks = peaks;
kymoseg_rec1.pair_ind_array = pair_ind_array;
kymoseg_rec1.separator = separator;
kymoseg_rec1.peak_linking_table = peak_linking_table;
kymoseg_rec1.param_array = param_array;
kymoseg_rec1.path = path;

cd(path.save);   
save('kymoseg_rec1', 'kymoseg_rec1', '-v7.3');



ret = 1;

end

% ===================================================================== %






