
% For segmenting the cell image by fluorescence channel
% Remove all small and boundary objects from the beginning.

function [output] = SegPhaseV8B(input, para)

% input = img0;
% para.conv_para = [3 1];          % (1) Convolution radius (2) Convolution sigma 
% para.adapt_para = [5 6000 0.2];  % (1) Adaptive threshold radius (2) minimal level (3) adaptive sigma 
% para.size_para = [80 1000];      % minimal/maximal cell size (in pixel)
% para.dilate_para = [3 3];        % Gaussian parameter for dilation (radius, SD in pixel)

% === (STEP 1) Segmentation the inner and outer contour ====== %

input = double(input);
img_conv = Img_GS_conv(input, para.conv_para);
mask1 = Mask_adap_thr(img_conv, para.adapt_para1);   % inner   
mask2 = Mask_adap_thr(img_conv, para.adapt_para2);   % outer

% === (STEP 2) filter small objects  ============= %

% inner contour (core region)
obj_seg1 = bwlabel(mask1); 
obj_seg1 = OJ_size_filter(obj_seg1, para.size_para);
obj_seg1 = OJ_reindex(obj_seg1);

output.mask1 = (obj_seg1 > 0);
output.mask1 = Mask_fill(output.mask1);

% outer contour (cell contour)
obj_seg2 = bwlabel(mask2); 
obj_seg2 = OJ_size_filter(obj_seg2, para.size_para);
obj_seg2 = OJ_reindex(obj_seg2);

output.mask2 = (obj_seg2 > 0);
output.mask2 = Mask_fill(output.mask2);

if (para.dilate_option == 1)

    output.mask2 = Img_GS_conv(output.mask2, para.dilate_para);
    output.mask2 = (output.mask2 > 0.1);

end

% === (STEP 3) Save images for print out  ========= %

colormap1 = gray;
colormap2 = gray;

OL_mask = double(output.mask1) + double(output.mask2);

output.Fig1 = RGB_format(img_conv, colormap1);
output.Fig2 = RGB_format(OL_mask, colormap2);

% figure;
% subplot(131); imagesc(img_conv);
% subplot(132); imagesc(OL_mask);
% axis equal;


