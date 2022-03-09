
% 2022 Feb 21th: Final version. This code is doing: 
% Segmentation for each debacked phase image.
% Two parameter sets are used: one for "core region", the other for "cell contour". 
% The "core region" is used for 2D kymograph segmentation in the next step,
% the cell outline is used to measure cell size and fluorescence intensity.

% Perform pre-segmentation on the debacked image. 
addpath(strcat(script_dir, 'B0_stack_preseg/segmentation_tools/'));

% Specified the subfolders
SubFolderPrefix = 'waf7_001xy3'; 
SubFolderList = [1:10]'; 

% ================================================================== %

paraA = {}; 
paraA.stack_range = [1 2000];   % Range of image stack    
paraA.file.prefix = strcat(SubFolderPrefix,'t');
paraA.file.postfix = 'c1DB.tif';
paraA.max_digit = max_digit;
paraA.flip_option = 0;

%--------------------------------------------------------------------% 
% ------------ Assign the segmentation parameter below ------------- %

% This code can set up different parameters for different frame.
% User must define the parameter for the 1st frame. For following frames:
% (1) If there is no new parameter, use the same parameter as the preious frame
% (2) If there is a new parameter, use the new parameter.

paraB = {};

% Create segmentation parameter for each frame
for fn = paraA.stack_range(1):paraA.stack_range(2)
    paraB{fn}.flag = 0;
end

m = 1;                                  % Frame number
paraB{m}.conv_para = [7 3];             % Gaussian convolution (radius, SD in pixel) for phase image
paraB{m}.adapt_para1 = [7 9000 0.2];    % Parameter for "core region"
paraB{m}.adapt_para2 = [7 8000 0.1];    % Parameter for "cell contour"
paraB{m}.size_para = [30 20000];        % Minimal/maximal cell size
paraB{m}.dilate_option = 1;             % Dilation flag
paraB{m}.dilate_para = [5 2];           % Dilate segmentation result by Gaussain convilution (radius, SD)
paraB{m}.flag = 2;                      % Flag for user-defined parameter

% m = 200;
% paraB{m}.conv_para = [5 3];
% paraB{m}.adapt_para1 = [7 10500 0.25];
% paraB{m}.adapt_para2 = [7 9700 0.1];
% paraB{m}.size_para = [30 20000];
% paraB{m}.dilate_option = 1; 
% paraB{m}.dilate_para = [5 2]; 
% paraB{m}.flag = 2;


% ------------ Assign the segmentation parameter above -------------- %
% =================================================================== %

parfor j = 1:size(SubFolderList,1)
    
    TagJ = num2str( SubFolderList(j));
    PathJ = strcat(data_folder, SubFolderPrefix, '-', TagJ, '\' );
    
    cd(script_dir);       
    [ret1] = sub_preseg(paraA, paraB, PathJ);
    
end

close all;


