function [cellMask] = phaseModule(A0,seg_param)

% current_folder = fileparts(mfilename('fullpath'));
% addpath(strcat(current_folder,'/segmentation_tools'));
% seg_param = [7 3500 0.2];
% (1) Load image
% A0 = double(imread('430_22may2018_pye_30c_001xy1c1.tif'));

A1 = Img_invert_shift(A0, 1, 5000);

wave_cut = 6;
A2 = Img_FFT_denoise(A1, wave_cut);

GS_param = [2 1];
A3 = Img_GS_conv(A2, GS_param);

B1 = Mask_adap_thr(A3, seg_param);

cellMask = bwlabel(B1);

figure;
subplot(131);imagesc(A3);
subplot(132);imagesc(B1);
subplot(133);imagesc(cellMask);