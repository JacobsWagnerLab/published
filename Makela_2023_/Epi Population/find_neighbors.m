%% Number of neighbours for each cell
% Calculates number of neighbors for a cell mask. This is done by first 
% dilating cell mask X number of pixels and seeing how many unique mask IDs
% overlap with the dilated mask.

% 4/27/2021 Jarno Makela

min_area = 150;             % min cell area
channel = 'c3.tif';         % ending of fluorescence channel e.g.'c3.tif'
dilate_size = 7;            % max neighbor distance from a cell in pixels

% find image mask files and extract XY number
images_mask_all = dir('*c1_mask.*');
for ii = 1:length(images_mask_all)
    newStr = split(images_mask_all(ii).name,["xy","c1_mask"]);
    images_mask_all(ii).XY = str2num(newStr{2});
end

cells = struct;         % this is where results are saved
cells.XY = [];              % XY field
cells.ID = [];              % mask ID
cells.no_neighb = [];       % number of neighbors
cells.neighb = [];          % mask IDs of neighbors
cells.area = [];            % area of a cell in pixels
cells.totalConc = [];       % concentration of a cell (per pixel)

n = 1;
uniq_XY = unique([images_mask_all.XY]);
for hh = 1:length(uniq_XY)  % loop over XY fields
    % image_mask files only for one XY
    images_mask = images_mask_all([images_mask_all.XY] == uniq_XY(hh));
    
    % load mask
    bw_mask = imread(images_mask.name);
    
    % ignore cells touching border
    bw_mask = imclearborder(bw_mask,4);

    % remove masks smaller than 150 pixels (4-connected) 
    % by multiplying the mask on the image
    bw_mask = bw_mask.*uint16(bwareaopen(bw_mask,min_area,4));

    % load fluorescence channel images
    % also filter fluorescence images with tophat to remove background
    newStr = split(images_mask.name,'c1_mask.tif');
    fl_I = imread([newStr{1} channel]);

    % cellIDs (remove 0 that is background)
    cellIDs = unique(bw_mask(:));
    cellIDs = double(cellIDs(cellIDs ~= 0));

    for ii = 1:length(cellIDs) % loop over masks
        % find mask area and dilate it by X pixels
        mask_cell = bw_mask == cellIDs(ii);
        mask_dilated = imdilate(mask_cell, strel('disk', dilate_size));
        
        % crop out area of dilated mask from image with all masks
        bw_mask_crop = double(bw_mask).*mask_dilated;
        % find out all unique mask IDs
        unique_masks = unique(bw_mask_crop);
        % remove self and 0 (background)
        unique_masks = unique_masks((unique_masks ~= 0) & (unique_masks ~= cellIDs(ii)));
        
        % store info to the struct
        cells(n).XY = uniq_XY(hh);
        cells(n).ID = cellIDs(ii);
        cells(n).no_neighb = length(unique_masks);
        cells(n).neighb = unique_masks;
        cells(n).area = sum(mask_cell(:));
        
        % calculate concentration by dividing sum of intensities by area
        fl_crop = double(fl_I).*mask_cell;
        cells(n).totalConc = sum(fl_crop(:))./sum(mask_cell(:));
        
        % update index
        n = n + 1;
    end
end
% vector with number of neighbors
vector = [cells.no_neighb];

% select specific number of neighbors
mask1 = vector == 0;
mask2 = vector == 1;
areas = [cells.area];
conc = [cells.totalConc];

figure;hold on;
plot(areas(mask1),conc(mask1),'.r')
plot(areas(mask2),conc(mask2),'.b')

[rho, pval] = corr(areas(mask1)',conc(mask1)','Type','Spearman')
