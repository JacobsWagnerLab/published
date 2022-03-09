
% Print substacks in one figure

function [ret1] = print_segmat_panel(segmatN, print_dir, param)

stack_range = param.Nrange;

image_per_panel = 60;
panel_num = ceil( (stack_range(2)-stack_range(1)+1)/ image_per_panel );

stack_index = {};

for k = 1:panel_num
    
    indA = stack_range(1) + (k-1) * image_per_panel ;
    indB = min(stack_range(2), indA + image_per_panel + 1  );  
    
    stack_index{k} = [indA indB];    
    
end

for k = 1:panel_num
    
    [ret1] = print_single_panel(segmatN, stack_index{k}, print_dir);

end

end

% ==================================================================== %

function [ret1] = print_single_panel(segmat, stack_range, print_dir)

original_seg = {};

rg = stack_range;

for j = rg(1):rg(2)
    
    original_seg{j-rg(1)+1} = segmat(:,:,j);

end

% Plot sub-stack data from SegMatContent

img_panel = {};

col_trim = 10; 
single_img_size = size(original_seg{1});
single_img_size(2) = single_img_size(2) - 2*col_trim;

img_num = length(original_seg); 

img_panel.seg_ori = NaN(single_img_size(1), single_img_size(2)*img_num);

for k = 1:img_num
    
    % Column range of single image for trimmimg
    col_trim_range = [col_trim col_trim + single_img_size(2)-1];
    
    % Column range on image panel
    col_range = (k-1)*single_img_size(2) + [1 single_img_size(2)];    

    % Paste each image on panel

    img_panel.seg_ori(:,col_range(1):col_range(2)) = ...
        original_seg{k}(:, col_trim_range(1):col_trim_range(2));       
    
end


%%%

index_ini = stack_range(1);

h1 = figure('Position', [1 1 size(img_panel.seg_ori, 2)  size(img_panel.seg_ori, 1)] );

imagesc('CData', img_panel.seg_ori);
set(gca,'LooseInset',get(gca,'TightInset'));

ax = gca;
ax.XAxisLocation = 'top';
ax.YDir = 'reverse';
axis equal;

xscale_points = single_img_size(2)* ( (1:5:img_num) -0.5);
xscale_index = (1:5:img_num) + index_ini - 1;
xticks(xscale_points);
xticklabels(xscale_index);

colormap(jet);

%

cd(print_dir);

pause(1);

save_name = strcat('seg_', num2str(stack_range(1)), '_', num2str(stack_range(2)),'.png' );
saveas(h1, save_name);

ret1 = 1;

close all;

end

% ==================================================================== %