
% Plot cell kymograph from a image stack

% ======================= Enter path variables ==========================

[filepath, ~, ~] = fileparts(mfilename('fullpath'));
addpath(strcat(filepath, '/segmentation_tools'));

main_path = '/Volumes/Data_04/Wei-Hsiang Lin/WHLin_data_N/(date)/(dataset)/';
subpath_c2 = 'c2Exp';
subpath_c3 = 'c3Exp';
filename_prefix = 'waf7_006xy1t';
filename_postfix_c2 = 'c2s.tif';
filename_postfix_c3 = 'c3s.tif';
subpathSeg = 'Seg/SegMatContent.mat';

% ================= Enter image analysis variables ======================

img_stack = {};

max_digit = 4;
Rsc_adjust_factor = 2;
index_range = [1 540];
frame_gap = 6;

% ======================= (1) Load images  ==============================

index_list = [];

for ind = index_range(1):index_range(2)
    if (rem(ind, frame_gap) == 1)
        index_list = [index_list; ind];
    end
end

segstack = load(strcat(main_path, subpathSeg)).SegMatContent;

for j = 1:length(index_list)
    
    index = index_list(j);
    
    filename_c2 = strcat(filename_prefix, AddZeros(max_digit,index), filename_postfix_c2);
    filename_c3 = strcat(filename_prefix, AddZeros(max_digit,index), filename_postfix_c3);
    
    filepath_c2 = strcat(main_path, subpath_c2, '/', filename_c2);
    filepath_c3 = strcat(main_path, subpath_c3, '/', filename_c3);
    
    img_stack{j}.index = index;
    img_stack{j}.c2 = imread(filepath_c2);
    img_stack{j}.c3 = imread(filepath_c3);        
    img_stack{j}.seg = segstack(:,:,index);
    
    img_stack{j}.Rsc_mat = image_ratio(img_stack{j}.c2, img_stack{j}.c3, img_stack{j}.seg, Rsc_adjust_factor);    
    
    
end

% ======================= (2) Concatenate images ==============================

Rsc_panel = [];
Rsc_trim_pixel = 15;

for j = 1:length(index_list)
    
    Rsc_temp = img_stack{j}.Rsc_mat;    
    Rsc_temp = Rsc_temp(:,1+Rsc_trim_pixel:end-Rsc_trim_pixel);
    Rsc_panel = [Rsc_panel Rsc_temp]; 
    
end

% ======================= (3) Plot the image ==============================

% Making customized RYB colormap with black background

RYBmap = customcolormap(linspace(0,1,11), {'#a60026','#d83023','#f66e44','#faac5d','#ffdf93','#ffffbd','#def4f9','#abd9e9','#73add2','#4873b5','#313691'}, 128);
RYB_gray = RYBmap;
RYB_gray(1,:) = 0.1*[1 1 1];

% Making plot

clim = [2 8];  
 
h1 = figure('position', [1 1 2520 421]);
imagesc(Rsc_panel, clim);
colormap(RYB_gray);

set(gca,'LooseInset',get(gca,'TightInset'));
colorbar('westoutside');
axis equal;
set(gca, 'YDir','normal');
set(gca,'xticklabel',[]);

single_img_size = size(Rsc_temp);
img_num = length(index_list);

xscale_points = single_img_size(2)* ( (0:5:img_num) - 0.5);
xscale_index = (0:5:img_num) + index_range(1)/frame_gap;
xticks(xscale_points);
xticklabels(frame_gap*xscale_index);

ylim([0 400]);

% ===================================================================== %

function [Rsc_mat] = image_ratio(imgA, imgB, seg, Rsc_adjust_factor)

mask = seg>0;
conv_param = [15 7];
A_conv = Img_GS_conv(imgA, conv_param);
B_conv = Img_GS_conv(imgB, conv_param);

bgA = mean(mean( imgA(1:10, 1:10) ));
bgB = mean(mean( imgB(1:10, 1:10) ));

% Trimming the mask (optional)

trim_level = 0.72;
trim_flag = 1;

if (trim_flag == 1)
    
    mask2 = trim_mask(mask, trim_level);

else
    
    mask2 = mask;

end

A_cell = (A_conv - bgA).*mask2;
B_cell = (B_conv - bgB).*mask2;

%

Rsc_mat = NaN(size(imgA,1), size(imgA,2));

for j1 = 1:size(imgA,1)
    
    for j2 = 1:size(imgA,2)
                        
        flag1 = ( A_cell(j1,j2) > 0 );
        flag2 = ( B_cell(j1,j2) > 0 );
        
        if ( flag1*flag2 == 1 )
            
            Rsc_mat(j1,j2) = A_cell(j1,j2)/B_cell(j1,j2);
            
        end
        
    end
    
end

Rsc_mat = Rsc_adjust_factor * Rsc_mat;


end

% ===================================================================== %

function [mask2] = trim_mask(mask, trim_level)

conv_param = [3 1];
mask2 = (Img_GS_conv(mask, conv_param) > trim_level);

end

% ===================================================================== %


function [StrOutput] = AddZeros(digit, j)
% Adding zero before number and generate string

% EXAMPLE
% digit = 3;
% j=132;

StrOutput = '';

% test how many digit for current string
StrTemp = num2str(j);
StrSize = size(StrTemp);    
ZeroCount = digit - StrSize(2);

if (ZeroCount>0)
    
    StrOutput = StrTemp;
    
    for m = 1:ZeroCount        
        StrOutput = strcat('0',StrOutput);        
    end
    
else
    
    StrOutput = StrTemp;
    
end

end

% ===================================================================== %