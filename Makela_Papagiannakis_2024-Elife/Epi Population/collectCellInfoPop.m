function collectCellInfoPop
% Collect cell information for a cells structure from population images.
% Cell areas are assumed to be saved into same folder as *c1_mask.* files,
% each mask having it's unique ID. Raw image data is included as cropped
% images for each cell.

%  Assumes the following filenames: 
%   ***xy**c*.tif
%       *c1.tif - phase contrast
%       *c1_mask.tif - mask file with each masking having unique number
%       *c2.tif - fluo1 channel
%       *c3.itf - fluo2 channel

% 5/4/21 Jarno Makela

no_fluo_channels = 2;           % number of fluorescence channels
size_min = 150;                 % minimum size for a cell in px   
pixel = 0.0666;                 % pixel size in um
padding = 6;                    % extra area (pixels in all directions) around cell for cropped images
dilate_size = 7;                % dilate mask to check neighbours along border
max_fraction = 0.05;            % maximum fraction of cells border occupied by other cells to be a considered a lone cell

% find image mask files and extract XY number
images_mask_all = dir('*c1_mask.*');
for ii = 1:length(images_mask_all)
    newStr = split(images_mask_all(ii).name,["xy","c1_mask"]);
    images_mask_all(ii).XY = str2num(newStr{2});
end

% initialize cells fields
cells = struct;
cells.XY = [];
cells.ID = [];
cells.coord = [];
cells.maskID = [];
cells.alone = [];
cells.neighb_IDs = [];
cells.mask = [];
cells.phase = [];
if no_fluo_channels == 1 
    cells.bg_fl1 = [];
    cells.fluo1 =[];
elseif no_fluo_channels == 2 
    cells.bg_fl1 = [];
    cells.bg_fl2 = [];
    cells.fluo1 = [];
    cells.fluo2 = [];
elseif no_fluo_channels == 3 
    cells.bg_fl1 = [];
    cells.bg_fl2 = [];
    cells.bg_fl3 = [];
    cells.fluo1 = [];
    cells.fluo2 = [];
    cells.fluo3 = [];
end

n = 1;
uniq_XY = unique([images_mask_all.XY]);
for hh = 1:length(uniq_XY)  % XY fields
    % image_mask files only for one XY
    images_mask = images_mask_all([images_mask_all.XY] == uniq_XY(hh));
    
    % load mask
    bw_mask = imread(images_mask.name);
    
    % background area
    bg_mask = bw_mask == 0;
    
    % ignore cells touching border
    bw_mask = imclearborder(bw_mask,4);

    % remove masks smaller than X (4-connected) 
    % by multiplying the mask on the image
    bw_mask = bw_mask.*uint16(bwareaopen(bw_mask,size_min,4));

    % load fluorescence channel images
    % estimate background intensity
    % also filter fluorescence images with tophat to remove background
    newStr = split(images_mask.name,'c1_mask.tif');
    if no_fluo_channels == 0
        ph_I = imread([newStr{1} 'c1.tif']);
    elseif no_fluo_channels == 1
        ph_I = imread([newStr{1} 'c1.tif']);
        fl1_I = imread([newStr{1} 'c2.tif']);
        bg_fl1 = median(fl1_I(bg_mask));
    elseif no_fluo_channels == 2
        ph = imread([newStr{1} 'c1.tif']);
        fl1_I = imread([newStr{1} 'c2.tif']);
        bg_fl1 = median(fl1_I(bg_mask));
        fl2_I = imread([newStr{1} 'c3.tif']);
        bg_fl2 = median(fl2_I(bg_mask));
	elseif no_fluo_channels == 3
        ph = imread([newStr{1} 'c1.tif']);
        fl1_I = imread([newStr{1} 'c2.tif']);
        bg_fl1 = median(fl1_I(bg_mask));
        fl2_I = imread([newStr{1} 'c3.tif']);
        bg_fl2 = median(fl2_I(bg_mask));
        fl3_I = imread([newStr{1} 'c4.tif']);
        bg_fl3 = median(fl3_I(bg_mask));
    end

    % cellIDs
    cellIDs = unique(bw_mask(:));
    cellIDs = double(cellIDs(cellIDs ~= 0));

    for ii = 1:length(cellIDs) % masks
        % cell mask from first time point
        mask_cell = bw_mask == cellIDs(ii);
        
        % find number of neighbours
        % dilate mask by X pixels
        mask_dilated = imdilate(mask_cell, strel('disk', dilate_size));
        % boundary of dilated mask
        boundary = bwboundaries(mask_dilated,8,'noholes');
        % mask values at the boundary
        pixels = impixel(bw_mask,boundary{1}(:,2),boundary{1}(:,1));
        % fraction of boundary having neighbors; use only first column
        fraction = sum(pixels(:,1) ~= 0)./length(pixels(:,1));
        % cell is considered to be alone based on the fraction
        if fraction > max_fraction
            cells(n).alone = 0;
            unique_masks = unique(pixels(:,1));
            unique_masks = unique_masks((unique_masks ~= 0));
            cells(n).neighb_IDs = unique_masks;
        else
            cells(n).alone = 1;
        end

        % crop image for each channel using padding around the mask
        % crop only if unique frames exists for channel
        s = regionprops(mask_cell,'BoundingBox');
        s = s(1);
        pd = padding;
        if no_fluo_channels == 0
            mask_crop = imcrop(mask_cell,s.BoundingBox+[-pd -pd 2*pd-1 2*pd-1]);
            ph_crop = imcrop(ph_I,s.BoundingBox+[-pd -pd 2*pd-1 2*pd-1]);
        elseif no_fluo_channels == 1
            mask_crop = imcrop(mask_cell,s.BoundingBox+[-pd -pd 2*pd-1 2*pd-1]);
            ph_crop = imcrop(ph_I,s.BoundingBox+[-pd -pd 2*pd-1 2*pd-1]);
            fl1_crop = imcrop(fl1_I,s.BoundingBox+[-pd -pd 2*pd-1 2*pd-1]);
        elseif no_fluo_channels == 2
            mask_crop = imcrop(mask_cell,s.BoundingBox+[-pd -pd 2*pd-1 2*pd-1]);
            ph_crop = imcrop(ph,s.BoundingBox+[-pd -pd 2*pd-1 2*pd-1]);
            fl1_crop = imcrop(fl1_I,s.BoundingBox+[-pd -pd 2*pd-1 2*pd-1]);
            fl2_crop = imcrop(fl2_I,s.BoundingBox+[-pd -pd 2*pd-1 2*pd-1]);
        elseif no_fluo_channels == 3
            mask_crop = imcrop(mask_cell,s.BoundingBox+[-pd -pd 2*pd-1 2*pd-1]);
            ph_crop = imcrop(ph,s.BoundingBox+[-pd -pd 2*pd-1 2*pd-1]);
            fl1_crop = imcrop(fl1_I,s.BoundingBox+[-pd -pd 2*pd-1 2*pd-1]);
            fl2_crop = imcrop(fl2_I,s.BoundingBox+[-pd -pd 2*pd-1 2*pd-1]);
            fl3_crop = imcrop(fl3_I,s.BoundingBox+[-pd -pd 2*pd-1 2*pd-1]);
        end

        % insert mask into cells
        cells(n).mask = mask_crop;
        cells(n).phase = ph_crop;
        
        % insert fluorescence channels
        if no_fluo_channels == 1
            cells(n).fluo1 = fl1_crop;
            cells(n).bg_fl1 = bg_fl1;
        elseif no_fluo_channels == 2 
            cells(n).fluo1 = fl1_crop;
            cells(n).fluo2 = fl2_crop;
            cells(n).bg_fl1 = bg_fl1;
            cells(n).bg_fl2 = bg_fl2;
        elseif no_fluo_channels == 3 
            cells(n).fluo1 = fl1_crop;
            cells(n).fluo2 = fl2_crop;
            cells(n).fluo3 = fl3_crop;
            cells(n).bg_fl1 = bg_fl1;
            cells(n).bg_fl2 = bg_fl2;
            cells(n).bg_fl3 = bg_fl3;
        end  
        
        % convert area from px into um and cells with less than min_frames
        cells(n).area = sum(mask_cell(:)).*pixel^2;
        cells(n).XY = uniq_XY(hh);
        cells(n).ID = n;
        cells(n).maskID = cellIDs(ii);

        % mask coordinates
        [y,x] = find(mask_cell);
        cells(n).coord = [mean(x) mean(y)];

        % update index
        n = n + 1;
    end
end

% save individual cells file
folderName = split(pwd,'\');
name_cells = char(strcat(folderName{end},'_cells.mat'));
save(name_cells,'cells')
end

