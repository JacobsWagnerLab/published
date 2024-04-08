% ======================================================================

% This scripts analyzes spot-markers such as Div-CFP and PleC-YFP data.
% The input is raw .tif images of phase-contrast, CFP and YFP channels. 

% The scripts call "double_marker" which is a pipeline to analyze makers
% and assign cell types based on the marker.

% ======================================================================

script_path = '(this_dir)';
addpath(strcat(script_path, 'segmentation_tools'));
addpath(strcat(script_path, 'cell_cord_tools'));

data_path = '(data_dir)';
data_subfolder = '(data_dir/subfolder)';
file_name = '(file_name)';
save_path = data_path;
save_subfolder = data_subfolder;

cd(strcat(data_path, data_subfolder));

%%% Read and save paths 

read_path = {};
read_path.c1 = strcat(data_path, data_subfolder, file_name,'c1.tif');
read_path.c2 = strcat(data_path, data_subfolder, file_name,'c2.tif');
read_path.c3 = strcat(data_path, data_subfolder, file_name,'c3.tif');

save_path = strcat(save_path, save_subfolder);

%%% Parameters for segmentation, spot detection and cell typing

seg_para = {};
seg_para.core_para = [7 5000 0];  %
seg_para.terra_param = [7 5];
seg_para.cell_param = [9 1000 0];
seg_para.mask_dilation = 1;
seg_para.size_param = [150 1000];

spot_para = {};
spot_param.c2 = [5 200 0];
spot_param.c3 = [5 200 0];

typing_param = {};

typing_param.c2.pole = 1;
typing_param.c2.pole_frac = 0.2;
typing_param.c2.mid = 0;
typing_param.c2.mid_frac = NaN;
typing_param.c2.minimal_size = 4;
typing_param.c2.minimal_avg_level = 200;

typing_param.c3.pole = 1;
typing_param.c3.pole_frac = 0.2;
typing_param.c3.mid = 0;
typing_param.c3.mid_frac = NaN;
typing_param.c3.minimal_size = 4;
typing_param.c3.minimal_avg_level = 300;

%%%

cd(script_path);
[retA] = double_markers(read_path, save_path, file_name, seg_para, spot_param, typing_param);


% ======================================================================
% ======================================================================

function [retA] = double_markers(read_path, save_path, file_name, seg_para, spot_param, typing_param)

% This scripts input raw images and perform cell typing based on
% two fluorescence markers. 


% (I-1) Read files 

img1 = double(imread(read_path.c1));
img2 = double(imread(read_path.c2));
img3 = double(imread(read_path.c3));

% (I-2) Use phase channel to do segmentation, and sort the segmented object
%       by area

save_choice = 0;
[seg_mat, BacData] = SegPhaseV6A(img1, seg_para, save_choice, '', '');
[seg_matS, BacDataS] = sub_sort_area(seg_mat, BacData);

% (I-3) Use segmentation mask to perform background subtraction

sampling_interval = 10;
[bg_function2, bg2] = ConstructBGFunction(img2, seg_matS, sampling_interval);
img2F = img2 - bg2;
img2F = max(img2F, 1);

sampling_interval = 10;
[bg_function3, bg3] = ConstructBGFunction(img3, seg_matS, sampling_interval);
img3F = img3 - bg3;
img3F = max(img3F, 1);

%  ====================================================================== %

% (II-1) Tile-display for all cells

m = size(BacDataS, 2);
crop_size = 45;
panel_col = 15;
panel_row = ceil(m/panel_col);

[cord_mat, tile_mat] = sub_get_coordinate(BacDataS, seg_matS, crop_size, panel_col);

tile = {};
tile.c1 = sub_crop_tile(img1, cord_mat, tile_mat, crop_size, panel_row, panel_col);
tile.c2 = sub_crop_tile(img2F, cord_mat, tile_mat, crop_size, panel_row, panel_col);
tile.c3 = sub_crop_tile(img3F, cord_mat, tile_mat, crop_size, panel_row, panel_col);
tile.seg = sub_crop_tile_filter(seg_matS, cord_mat, tile_mat, crop_size, panel_row, panel_col);

tile.seg_para = seg_para;
tile.cord_mat = cord_mat;
tile.tile_mat = tile_mat;

% (II-2) For each object on Tiles image, perform cell-coordinate analysis

tile.seg_data = regionprops(tile.seg, 'Centroid','Orientation','Area');
tile.cell_cord = cell_coordinate(tile.seg);

%  ====================================================================== %

% (III-1) Spot analysis on background-subtracted and tiled images.

spot_data2 = {};   
spot_data2.param = spot_param.c2;
[spot_data2.ensemble] = spot_analysis(tile.c2, spot_data2.param);

spot_data3 = {};   
spot_data3.param = spot_param.c3;
[spot_data3.ensemble] = spot_analysis(tile.c3, spot_data3.param);

% (III-2) Associate each cell with closed spots

[spot_data2.cell_cord, spot_data2.spot2cell] = associate_cell_to_spot(tile, spot_data2.ensemble);

[spot_data3.cell_cord, spot_data3.spot2cell] = associate_cell_to_spot(tile, spot_data3.ensemble);

%  ====================================================================== %

% (IV) Cell-typing by spot criteria

[cell_type] = cell_typing_GRIT2(tile, spot_data2, spot_data3, param);

%%%  ================================================================= %

% (V-1) Saving images

cd(save_path);
[ret1] = sub_save_tile_image_M2(tile, spot_data2, spot_data3, cell_type, file_name);

% (V-2) Saving data

cd(save_path);

data = {};
data.tile = tile;
data.spot_data2 = spot_data2;
data.spot_data3 = spot_data3;
data.cell_type = cell_type;

save(strcat(file_name, '_data.mat'), 'data');
retA = 1;


% ======================================================================
% ======================================================================

% Subroutines that used in the double-marker pipeline:
% - sub_sort_area
% - sub_get_coordinate
% - sub_crop_tile
% - sub_crop-tile_filter
% - spot_analysis
% - associate_cell_to_spot
% - cell_cyping_GRIT2
% - sub_save_tile_image_M2

% ======================================================================
% ======================================================================

% Sort BacDataF by object area

function [seg_matS, BacDataS] = sub_sort_area(seg_mat, BacData)

m = size(BacData,1);

area = NaN(m,2);

for j = 1:m
    
    area(j,:) = [BacData(j).Area j];

end

area_sorted = sortrows(area,1);

for j = 1:m
   
   ind = area_sorted(j,2);
   BacDataS(j) = BacData(ind);

end 

% Reploting segmentation matrix

seg_matS = 0*seg_mat;

% finding the reverse index
index_new_old = [area_sorted(:,2) (1:m)'];
index_new_oldS = sortrows(index_new_old,1);

for r = 1:size(seg_mat,1)
    for c = 1:size(seg_mat,2)
        
        if (seg_mat(r,c) > 0)  
            
            ind = index_new_oldS(seg_mat(r,c), 2);
            seg_matS(r,c) = ind;

        end
        
    end
end

% ======================================================================
% ======================================================================

% Get coordinate for cell tiling

function [cord_mat, tile_mat] = sub_get_coordinate(data, seg_mat, crop_size, panel_col)

m = size(data, 2);
mx = size(seg_mat);

centroid_mat = NaN(m,2);
cord_mat = NaN(m,4);
tile_mat = NaN(m,4);

for j = 1:m
    
    centroid_mat(j,:) = round(data(j).Centroid);
    
    % Note that {x,y} = {c,r} and the image is flipped
    
    cord_mat(j,1) = centroid_mat(j,2) - crop_size;
    cord_mat(j,2) = centroid_mat(j,2) + crop_size - 1;
    cord_mat(j,3) = centroid_mat(j,1) - crop_size;
    cord_mat(j,4) = centroid_mat(j,1) + crop_size - 1;

    cord_mat(j,1) = max(cord_mat(j,1), 1);
    cord_mat(j,2) = min(cord_mat(j,2), mx(2));
    cord_mat(j,3) = max(cord_mat(j,3), 1);
    cord_mat(j,4) = min(cord_mat(j,4), mx(1));
    
end

unit_L = 2*crop_size;

% When making the tile, filled up all space in the 1st row first,
% and move to the next row.

for j = 1:m
    
    tile_row_num = ceil(j/panel_col);
    tile_col_num = j - (panel_col * (tile_row_num-1));
    
    %[tile_row_num tile_col_num]
    
    tile_mat(j,1) = (tile_row_num - 1) *unit_L + 1;
    tile_mat(j,2) = tile_mat(j,1) + (cord_mat(j,2)-cord_mat(j,1)+1) - 1; 
    
    tile_mat(j,3) = (tile_col_num - 1) *unit_L + 1;
    tile_mat(j,4) = tile_mat(j,3) + (cord_mat(j,4)-cord_mat(j,3)+1) - 1; 
      
end


% ======================================================================
% ======================================================================

% Crop tiles based on coordinates

function [tile_img] = sub_crop_tile(ori_img, cord_mat, tile_mat, crop_size, panel_row, panel_col)

unit_L = 2*crop_size;
tile_img = zeros(unit_L*panel_row, unit_L*panel_col);

for j = 1:size(cord_mat,1)
    
    % Obtain the cropping region
    
    region_crop = ori_img(cord_mat(j,1):cord_mat(j,2), ...
                          cord_mat(j,3):cord_mat(j,4));
                      
    tile_img(tile_mat(j,1):tile_mat(j,2), ...
             tile_mat(j,3):tile_mat(j,4)) = region_crop;
        
end

% ======================================================================
% ======================================================================

% Crop tiles based on coordinates, filter with object index

function [tile_img] = sub_crop_tile_filter(ori_img, cord_mat, tile_mat, crop_size, panel_row, panel_col)

unit_L = 2*crop_size;
tile_img = zeros(unit_L*panel_row, unit_L*panel_col);

for j = 1:size(cord_mat,1)
    
    % Obtain the cropping region
    
    region_crop = ori_img(cord_mat(j,1):cord_mat(j,2), ...
                          cord_mat(j,3):cord_mat(j,4));
                      
    region_crop = region_crop .* (region_crop == j);
                      
    tile_img(tile_mat(j,1):tile_mat(j,2), ...
             tile_mat(j,3):tile_mat(j,4)) = region_crop;
        
end

% ======================================================================
% ======================================================================

% Characterize spot and perform spot statistics

function spot_ensemble = spot_analysis(img, spot_param)

%img = tile.c2;
%spot_param = [5 1000 0];

% Convolve with a narrow Gaussian
A1 = img;

GS_param = [3 1];
A2 = Img_GS_conv(A1, GS_param);

% Adaptive threshold
A3 = Mask_adap_thr(A2, spot_param);
spot_label = bwlabel(A3);
spot_data = regionprops(spot_label, 'Centroid','Orientation','Area');

% Analyze area, total inttensities of each spot

spot_num = size(spot_data,1);
spot_statistics = ones(spot_num, 3);  % area, total intensity, mean intensity

for r = 1:size(spot_label,1)
    
    for c = 1:size(spot_label,2)
        
        index = spot_label(r,c);
        
        if ( index > 0 )
            
            spot_statistics(index,1) = spot_statistics(index,1) + 1;
            spot_statistics(index,2) = spot_statistics(index,2) + A1(r,c);
            
        end
        
    end
    
end

spot_statistics(:,3) = spot_statistics(:,2) ./ spot_statistics(:,1);

%%%

spot_ensemble = {};
spot_ensemble.avg_size = mean(spot_statistics(:,1));
spot_ensemble.std_size = std(spot_statistics(:,1));
spot_ensemble.avg_meanI = mean(spot_statistics(:,3));
spot_ensemble.std_meanI = std(spot_statistics(:,3));

spot_ensemble.stat = spot_statistics;
spot_ensemble.label = spot_label;
spot_ensemble.location = spot_data;

% ======================================================================
% ======================================================================

% Associate spot with cells

function [record_cell, record_spot2cell] = associate_cell_to_spot(tile, spot_ensemble)

cell_seg = tile.seg;             % index-labeled images
cell_data = tile.seg_data;       % cell object data

spot_seg = spot_ensemble.label;  % index-labeled images
spot_location = spot_ensemble.location;   % spot object data

%

cell_num = max(size(cell_data));
spot_num = max(size(spot_location));

max_dist = 50;
record_spot2cell = zeros(spot_num, 1);

% Note that here, cell and spot location are in Centroid coordinate.
% {x,y} = {c,r}

for c = 1:cell_num
    
    cell_loc = cell_data(c).Centroid;
    
    for p = 1:spot_num
       
        spot_loc = spot_location(p).Centroid;
        
        dist_temp = sqrt(sum(cell_loc - spot_loc).^2);
        
        if ( dist_temp < max_dist )
            
            % Check for overlap betwen spot and "dilated cell".
            % If there is at least one pixel OL, 
            % associate this spot p with cell c.
            
            % Note that here we assume one spot associate to only one cell.
            % But one cell can be associated to multiple spot
            
            cell_range = (cell_seg == c);
            dilated_range = ( conv2(cell_range, PM_2D_gaussian(5,3), 'same') );
            dilated_range = (dilated_range > 0);  % binary mask
            
            spot_range = (spot_seg == p);
            
            OL_pixel = sum(sum(dilated_range.*spot_range));
            
            if (OL_pixel > 1)
                
                record_spot2cell(p) = c;
                
            end
            
        end
                       
        
    end
    
end
    
% For each cell, find all spot associated with the cell.

record_cell = {};

for c = 1:cell_num
    
    flag = 1;
    spot_collect = {};
    
    for j = 1:spot_num
       
        if ( record_spot2cell(j) == c)            
            
            spot_collect{flag}.spot_index = j;
            [ind_flag_percentile] = map_spot_on_cell_cord(c, j, tile, spot_ensemble);
            
            spot_collect{flag}.spot_midline_percentile = ind_flag_percentile; 
            
        end
        
    end
    
    record_cell{c} = spot_collect;
    
end

% ======================================================================
% ======================================================================

% Cell typing

function [cell_type] = cell_typing_GRIT2(tile, spot_data2, spot_data3, param)

% Cell typing. 
% This code based on a collection of criteria on spot property and cell
% size to classify the cell type
% param = typing_param;

cell_num = size(tile.tile_mat,1);

check2 = {};

check2.spot_c2_present = 1;
check2.spot_c2_min_size = 1;
check2.spot_c2_pole = 1;

%

check3.spot_c3_present = 1;
check3.spot_c3_min_size = 1;

% ============================================================= %

% Checking c2 spot

c2_criteria = zeros(cell_num, 1);

for j = 1:cell_num
    
    spot_dataJ = spot_data2.cell_cord{j};   
    spot_num = max(size(spot_dataJ)); 
    
    if ( spot_num == 0 )    % (1) Test if spot presence
        
        c2_criteria(j) = 0;
        
        
    elseif ( spot_num > 0 )  % (2) Test for each spot 
               
        spot_typing_matrix = NaN(spot_num, 2);  % Testing two things: pole, intensity
            
        for k = 1:spot_num  % check each spot
              
            % (i) Test if spot is at pole
            
            temp = spot_dataJ{k}.spot_midline_percentile;
            pole_test = ( temp < param.pole_frac(2) ) || ( temp > (1 - param.pole_frac(2)) ); 
            spot_typing_matrix(k,1) = pole_test;
            
            % (ii) Test if spot exceed minimal size
                        
            spot_focus = spot_dataJ{k};      
            spot_focus_size = spot_data2.ensemble.stat(spot_focus.spot_index, 1);    
            spot_typing_matrix(k,2) = (spot_focus_size > param.spot_c2_minimal_size);              
            
         end    
    
        % (iii) Checking if there is a spot satisfies all criteria
    
        checking_flag = 0;
    
        for k = 1:spot_num        
            checking_flag = checking_flag + (sum(spot_typing_matrix(k,:)) == 2);         
        end
    
        c2_criteria(j) = (checking_flag >=1 );
    
    end
    
end


% =================================================================== %

% Checking c3 spot

c3_criteria = zeros(cell_num, 1);

for j = 1:cell_num
    
    spot_dataJ = spot_data3.cell_cord{j};   
    spot_num = max(size(spot_dataJ)); 
    
    if ( spot_num == 0 )    % (1) Test if spot presence
        
        c3_criteria(j) = 0;
        
        
    elseif ( spot_num > 0 )  % (2) Test for each spot 
               
        spot_typing_matrix = NaN(spot_num, 1);  % Testing spot intensity
            
        for k = 1:spot_num  % check each spot
                          
            % (i) Test if spot exceed minimal size
                        
            spot_focus = spot_dataJ{k};      
            spot_focus_size = spot_data3.ensemble.stat(spot_focus.spot_index, 1);    
            spot_typing_matrix(k,1) = (spot_focus_size > param.spot_c3_minimal_size);              
            
         end    
    
        % (ii) Checking if there are TWO spots satisfying all criteria
    
        checking_flag = 0;
    
        for k = 1:spot_num        
            checking_flag = checking_flag + (spot_typing_matrix(k,:) == 1);         
        end
    
        c3_criteria(j) = (checking_flag >=2 );
    
    end
    
end


% ============================================================= %

% Cell typing

cell_type = zeros(cell_num,1);


for j = 1:cell_num
    
    PN_flag = [c2_criteria(j) c3_criteria(j)];
    
    if ( PN_flag == [1 0] ) 
        
        cell_type(j) = 1;
        
    elseif ( PN_flag == [0 0] ) 
        
        cell_type(j) = 2; 
            
    elseif ( PN_flag == [0 1] ) 
        
        cell_type(j) = 3; 
        
    elseif ( PN_flag == [1 1] ) 
        
        cell_type(j) = 4; 

    end
    
end

% ======================================================================
% ======================================================================

% Save tile images

function [ret1] = sub_save_tile_image_M2(tile, spot_data2, spot_data3, cell_type, file_name)

% ================= (1) Save the tiled raw images ========================

multiplier = [0.5 1.5 1];

OLmat = zeros(size(tile.c1,1), size(tile.c1,2), 3);  % R,G,B

c1 = multiplier(1) * tile.c1 / max(max(tile.c1));
c2 = multiplier(2) * tile.c2 / max(max(tile.c2));
c3 = multiplier(3) * tile.c3 / max(max(tile.c3));

% Set c1 to gray;
OLmat(:,:,1) = OLmat(:,:,1) + c1;
OLmat(:,:,2) = OLmat(:,:,2) + c1;
OLmat(:,:,3) = OLmat(:,:,3) + c1;

% Set c2 to red
OLmat(:,:,1) = OLmat(:,:,1) + c2;
%OLmat(:,:,2) = OLmat(:,:,2) + c2;
%OLmat(:,:,3) = OLmat(:,:,3) + c2;

% Set c3 to green
%OLmat(:,:,1) = OLmat(:,:,1) + c3;
OLmat(:,:,2) = OLmat(:,:,2) + c3;
%OLmat(:,:,3) = OLmat(:,:,3) + c2;

%figure; image(OLmat); axis equal; axis off;
imwrite(OLmat, strcat(file_name, '_OL.png'), 'png');

% ================= (2) Save the tiled segmented image ==================== 

h1 = figure; %('position', [1 1 size(tile.c1,2) size(tile.c1,1)]);
set(gca,'LooseInset',get(gca,'TightInset'));

spot_representation ...
    = 60*(spot_data2.ensemble.label > 0) + 30*(spot_data3.ensemble.label > 0) + 5*(tile.seg > 0);

image(spot_representation); hold on;

cm_black = jet;
cm_black(1,:) = 0.5*[1 1 1];
colormap(cm_black);  %axis equal; 
%axis equal; 
axis off;

tile_mat = tile.tile_mat;

for r = 1:size(cell_type,1)
    
    c23_string = '';
    
    if (cell_type(r) == 1)
        c23_string = '+,-';
    elseif (cell_type(r) == 2)
        c23_string = '-.+';
    elseif (cell_type(r) == 3)
        c23_string = '+.+';
    elseif (cell_type(r) == 4)
        c23_string = '-.-';
    end  
    
    type_string = strcat(c23_string);
    
    ind1 = tile_mat(r,1);
    ind2 = tile_mat(r,3);
    text(ind2, ind1 + 5, type_string, 'fontsize', 10, 'color', [1 1 1]); hold on;

end

hold off;

output_size = [size(tile.c1,2) size(tile.c1,1)];  %Size in pixels
resolution = 200;  %Resolution in DPI
set(h1,'paperunits','inches','paperposition',[0 0 output_size/resolution]);

tile_save_name = strcat(file_name, '_Typing.png');
print(tile_save_name,'-dpng',['-r' num2str(resolution)]);


% ================= (3)  Load two .png files and combine them (vertically) into a single image

OL_image = imread(strcat(file_name, '_OL.png'));
Typing_image = imread(strcat(file_name, '_Typing','.png'));

combined_image = zeros(2*size(OL_image,1), size(OL_image,2), 3);
rL = size(OL_image,1);
combined_image(1:rL, :,:) = OL_image;
combined_image(rL+1:2*rL, :,:) = Typing_image;

imwrite(uint8(combined_image), strcat(file_name, '_combined.png'), 'png');

% ================= (4) Delete the files of OL image and Typing image =====

delete(strcat(file_name, '_OL.png'));
delete(strcat(file_name, '_Typing.png'));

close all;
ret1 = 1;

% ======================================================================
% ======================================================================

