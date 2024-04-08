
% Given a segmentation matrix, get bounding box and cell coordinate for each objects

%function [obj_data] = cell_coordinate(seg_mat)
seg_mat = tile.seg;

obj_num = max(max(seg_mat));
obj_data = {};

m =85;
    
obj = {};
obj_location = (seg_mat == m);

temp_struct = regionprops(obj_location, 'BoundingBox');

% Note: switching {y,x} to {r,c}
box = temp_struct.BoundingBox;
obj.corner = floor([box(2) box(1)]');
obj.dim = ceil([box(4) box(3)]');

obj_pixel = obj_location(obj.corner(1): obj.corner(1)+obj.dim(1)-1,...
                         obj.corner(2): obj.corner(2)+obj.dim(2)-1);


% (1) assign level to obj_pixel
max_level = 15;
[obj.level] = find_level_mat(obj_pixel, max_level);

% (2) finding the midline
[obj.midline, mid_axis_output] = find_middle_axis(obj.level);

obj_data{m} = obj;



figure; 
imagesc(obj.level);  hold on;
plot(mid_axis_output(:,2), mid_axis_output(:,1), '-bo'); hold on;
plot(obj.midline(:,2), obj.midline(:,1), 'r-'); hold off;
colormap(gray);

%figure;  imagesc(obj_pixel);

%obj.centroid = temp_struct.Centroid;
%obj.area = temp_struct.Area;
%}
