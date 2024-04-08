
function [output, BacDataF] = SegPhaseV6A(input, seg_para, save_choice, save_folder, numtag)

% This protocol use phase image to segment cells,
% with 2-step algorithm
% 
% input = imread('test_imageDB.tif');
% A0 = input;
% save_folder = pwd;
% numtag = 1020;
% seg_para = {};

%seg_para.core_para = [7 10000 0.2];
%seg_para.terra_param = [7 5];
%seg_para.cell_param = [7 8000 0.1];
%seg_para.size_param = [80 1000];


% ================================================ %
% =====  (STEP 1) smoothing image           ====== %
% ================================================ %

A0 = Img_invert_shift(input, 1, seg_para.invert_scale);
A0 = max(1, A0);

% noise filter and smoothing image
wave_cut = 6;
A1 = Img_FFT_denoise(A0, wave_cut); 

GS_param = [3 1];
A2 = Img_GS_conv(A1, GS_param);

% ================================================ %
% =====  (STEP 2) generate terratery mask   ====== %
% ================================================ %

core_param = seg_para.core_para; 
terra_param = seg_para.terra_param;

B1 = Mask_adap_thr(A2, core_param);
B2 = Mask_terra(B1, terra_param);

% =========================================== %
% =====  (STEP 3) generate cell mask    ===== %
% =========================================== %

cell_param = seg_para.cell_param;

C1 = Mask_adap_thr(A2, cell_param); 
C2 = C1.*B2;

% ================================================= %
% =====  (STEP 4) filter bd and small objects  ==== %
% ================================================= %

size_param = seg_para.size_param;
bd_param = [2 3];

D1 = OJ_size_filter(C2, size_param);
D2 = OJ_bd_filter(D1, bd_param);
D3 = OJ_reindex(D2);

output = D3;
BacDataF = regionprops(output, 'Centroid','Orientation','Area');

% ================================================= %
% =====  (STEP 5) save png images  =========== ==== %
% ================================================= %

%figure; imagesc(B1);
%figure; imagesc(D3);

if (save_choice == 1)

colormap1 = hot;
colormap2 = colorcube;

Fig1 = RGB_format(input, colormap1);
Fig2 = RGB_format(B2, colormap2);
Fig3 = RGB_format(C2, colormap2);    
Fig4 = RGB_format(D3, colormap2);

img_size = size(input);
WImg = zeros(img_size(1), 4*img_size(2), 3);
    
for c = 1:3
    WImg(:,:,c) = [Fig1(:,:,c) Fig2(:,:,c) Fig3(:,:,c) Fig4(:,:,c)];
end

%%%

cd(save_folder);
WImgName = strcat('Figure', num2str(numtag),'.png');
imwrite(WImg, WImgName, 'png');    
close all;

end




