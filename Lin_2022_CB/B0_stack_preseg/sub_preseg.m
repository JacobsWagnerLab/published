
function [ret1] = sub_preseg(paraA, paraB, PathJ)

% (0) Get parameters
max_digit = paraA.max_digit;
stack_range = paraA.stack_range;
file = paraA.file;

% (1) Fill in parameter for segmentation on every frame
paraB = complete_seg_para(paraB, stack_range);


% (2) Perform pre-segmentation for each frame

cd(strcat(PathJ, '/c1Exp/'));
mask_data = {};

for m = stack_range(1):stack_range(2)
    
    filename = strcat(file.prefix, AddZeros(max_digit, m), file.postfix);
    img = imread(filename);
    
    if (paraA.flip_option == 1)  % perform flip on each column (i.e. y-axis)   
        img = flip(img);
    end
    
    mask_data{m}= SegPhaseV8B(img, paraB{m});    
    
end

% (3) Save results: 

cd(PathJ);
mkdir('PreSegB'); save_path = strcat(PathJ, 'PreSegB');

% (3-1) save figure files

[img_rgb] = save_figure_panel(mask_data, save_path, stack_range);

% (3-2) save data and parameter files

mask_record = {};

image_size = size(mask_data{1}.mask1);
stack_size = stack_range(2) - stack_range(1) + 1;

mask_record.SegMatContent = NaN(image_size(1), image_size(2), stack_size);
mask_record.SegMatContentD = NaN(image_size(1), image_size(2), stack_size);
mask_record.SegMatInd = NaN(stack_size);

for m = stack_range(1):stack_range(2)
    
    ind = m - stack_range(1) + 1;
    
    mask_record.SegMatContent(:,:,ind) = mask_data{m}.mask1;
    mask_record.SegMatContentD(:,:,ind) = mask_data{m}.mask2;
    mask_record.SegMatInd = m;
   
end

mask_record.paraA = paraA;
mask_record.paraB = paraB;
mask_record.img_rgb = img_rgb;

cd(save_path);
save_name = strcat('mask_record.mat');
save(save_name, 'mask_record', '-v7.3');

ret1 = 1;

end


% ======================================================================

function [img_rgb] = save_figure_panel(mask_data, save_path, stack_range)

trim_width = 15;
panel_max_frame = 30;

index_list = [];
index_gap = 15;

% (i) Determine how many image panels are displayed

for m = stack_range(1):stack_range(2)
    
    if ( rem(m, index_gap) == 1 )
        
        index_list = [index_list; m];
    end
    
end

index_list_array = {};

panel_num = ceil(length(index_list)/panel_max_frame);

for k = 1:panel_num
    
    temp_ind(1) = (k-1)*panel_max_frame + 1;
    temp_ind(2) = min(k * panel_max_frame, length(index_list));    
    index_list_array{k} = index_list(temp_ind(1):temp_ind(2));
    
end


% (ii) Concatenate the [debacked image] and [presegmentation data]

img_rgb = {};   % Image stacks. Each one is one panel

for k = 1:panel_num

index_vec = index_list_array{k};
    
img_mat = {};    
img_mat{1} = [];    
img_mat{2} = [];    
img_mat{3} = [];
    
for r = 1:length(index_vec)

    for c = 1:3
        
        temp1 = [];   % Debacked raw data
        temp2 = [];   % Pre-segmented data
        
        temp1 = mask_data{index_vec(r)}.Fig1(:,:,c);
        temp2 = mask_data{index_vec(r)}.Fig2(:,:,c);
    
        temp1 = temp1(:, trim_width : end-trim_width+1);
        temp2 = temp2(:, trim_width : end-trim_width+1);
        
        img_mat{c} = [img_mat{c} temp1 temp2];
    
    end
    
end

img_temp = NaN(size(img_mat{1},1), size(img_mat{1},2), 3 );
    
for c = 1:3
    img_temp(:,:,c) = img_mat{c};
end

img_rgb{k} = img_temp;

end


for k = 1:panel_num
    
h1 = figure('position', [1 1 size(img_rgb{k},2) size(img_rgb{k},1)]); 

% Plot image panel
imagesc(img_rgb{k});

colormap(gray);

axis equal;

unit_size = 2* size(temp1, 2);
unit_num = size(index_list_array{k}, 1);

xt = unit_size* ((1:unit_num)-0.5) + 1 ;
xtlbl = index_list_array{k};

set(gca, 'XTick',xt, 'XTickLabel',xtlbl);
set(gca,'LooseInset',get(gca,'TightInset'));


cd(save_path);

numA = index_list_array{k}(1);
numB = index_list_array{k}(end);
save_name = strcat('segmat_', num2str(numA), '-', num2str(numB),'.tif');

saveas(h1, save_name);

%close all;

end

ret1 = 1;

end

% ======================================================================

function para = complete_seg_para(para, stack_range)

% Use the user-input data, fill the segmentation parameter for each frame

for fn = ( stack_range(1) + 1) :stack_range(2)
    
    % Checkt each frame if new segmentation prameter is specified (para{}.flag == 2)
    % If not (para{}.flag == 0), use the same parameter as the previous frame
    % and mark para{}.flag = 1.

    if (para{fn}.flag == 0) && (para{fn-1}.flag > 0)
        para{fn} = para{fn-1};
        para{fn}.flag = 1;
    end
    
end

end

% ======================================================================

function [StrOutput] = AddZeros(digit, j)

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

% ======================================================================

