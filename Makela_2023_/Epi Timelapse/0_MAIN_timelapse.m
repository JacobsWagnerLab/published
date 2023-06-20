%% Time-lapse analysis
% Detect cell areas using U-net (or some other algorithm)
% In each frame, each cell mask should have a unique mask ID.

% Assumes the following filename structure: 
%  ***t**xy**c*.tif
%       *c1.tif - phase contrast
%       *c1_mask.tif - mask file with unique mask numbers
%       *c2.tif - fluo1 channel
%       *c3.itf - fluo2 channel
%       ...

% 7/7/22 Jarno Makela

%% Track cells and collect cell based image data
% track cells over consecutive frames using mask overlap and gives them
% unique cell IDs. Saves info on cell lineages to cell_info.mat
trackCells

% collect all image data into a .mat file according to cell tracking
% See options of the function
collectCellInfoTS

%% Tool to check correct segmentation
% assign 1 to correctly segmented or 0 to incorrectly segmented

min_frames = 10;                % minimum number of frames per cell
min_fold_size = 1.1;            % minimum area growth in fold change
area_last = 6;                  % minimum area in um of the last frame

% get cells file
[filename,pathname] = uigetfile('*.mat', 'Select','Select','MultiSelect', 'on');
if ~iscell(filename)
    name = filename;
    filename = cell(1);
    filename{1} = name;
end

figure('units','normalized','outerposition',[0 0 1 1])
for hh = 1:length(filename)
    load([pathname filename{hh}])

    cells(1).no_chrom = [];
    for ii = 1:length(cells)
        % enough frames and no no_chrom yet
        if length(cells(ii).area) > min_frames && isempty(cells(ii).no_chrom)
            if cells(ii).area(end)/cells(ii).area(1) > min_fold_size && ...
                    cells(ii).area(end) > area_last

                % choose frames - 4 frames distributed
                if isfield(cells,'fluo1')
                    fluo_exists = find(cellfun('isempty',cells(ii).fluo1) == 0);
                    indeces = ceil(linspace(1,length(fluo_exists),4));
                    frames = fluo_exists(indeces);
                else
                    frames = ceil(linspace(1,length(cells(ii).phase),4));
                end

                t = tiledlayout(2,4);
                % loop to plot images
                
                for jj = 1:length(frames)
                    % load images
                    bw_mask = cells(ii).mask{frames(jj)};
                    phase_image = cells(ii).phase{frames(jj)};
                    if isfield(cells,'fluo1')
                        fluo1_image = cells(ii).fluo1{frames(jj)};
                    else
                        fluo1_image = [];
                    end

                    % plotting
                    % phase contrast with auto gain
                    nexttile(jj);
                    imshow(imadjust(phase_image))
                    if ~isfield(cells,'fluo1')
                        hold on;
                        bw_mask_dil = imdilate(bw_mask, strel('disk', 2));
                        visboundaries(bw_mask_dil,'Color','r','LineWidth',0.25);
                    end

                    % HU fluoresence with auto gain
                    nexttile(jj+4);
                    if isfield(cells,'fluo1')
                        imshow(imadjust(fluo1_image))
                        hold on;
                        bw_mask_dil = imdilate(bw_mask, strel('disk', 2));
                        visboundaries(bw_mask_dil,'Color','r','LineWidth',0.25);
                    else
                        imshow([])
                    end
                    
                end
                % link axes and adjust padding
                t.TileSpacing = 'none';
                t.Padding = 'none';

                % input for number of chromosomes
                disp(['XY' num2str(cells(ii).XY)])
                disp(num2str(ii))
                str = input('Number of chromosomes or not suitable (0): ','s');  
                str = str2num(str);
                if ~isempty(str)
                    cells(ii).no_chrom = str; 
                else
                    cells(ii).no_chrom = 0;
                end
            else
                cells(ii).no_chrom = 0;
            end
        end
    end
    
    % save with same filename
    save([pathname filename{hh}],'cells','-v7.3')
end

% Remove unused cells from data files and join different XYs
joinCells

%% Truncate single cell time series
% truncates individual cell traces if instantaneous growth reduces below 
% maximum growth rate.

min_frames = 10;            % minimum number of frames per cell
win_size = 7;               % window size for smoothing growth rate
growth_fold_change = 0.75;  % cut threshold from maximum growth rate
min_streak = 3;             % how many frames growth has to be below threshold

% get cells file
[filename,pathname] = uigetfile('*.mat', 'Select','Select','MultiSelect', 'on');
if ~iscell(filename)
    name = filename;
    filename = cell(1);
    filename{1} = name;
end

for hh = 1:length(filename)
    load([pathname filename{hh}])

    cells(1).area_crop_est = [];
    for jj = 1:length(cells)
        % only cells with enough frames
        if length(cells(jj).area) > min_frames && ~isempty(cells(jj).no_chrom)     
            vector = cells(jj).area;
            vector = diff(vector);
            % smoothen the growth rate vector
            smooth_vector = movmean(vector,win_size);

            % find regions where before is growth_fold_change of max
            [~,max_frame] = max(smooth_vector);
            mask = smooth_vector./max(smooth_vector) < growth_fold_change;
            % determine frames before max not to be cropped
            mask(1:(max_frame)) = 0;

            % remove growth_fold_change streaks shorter than min_streak
            X = diff(mask~=0)~=0;
            B = find([true,X]); % begin of each group
            E = find([X,true]); % end of each group
            matrix = [];
            matrix(1,:) = B;
            matrix(2,:) = 1+E-B; % the length of each group
            matrix(3,:) = mask(B); % indeces of each group     
            matrix = matrix(:,matrix(3,:)==1); % only streaks of 1
            matrix = matrix(:,matrix(2,:)>=min_streak); % long enough
            
            % pick first 1 streak with enough indeces
            if ~isempty(matrix)
                cells(jj).area_crop_est = cells(jj).area(1:matrix(1,1));
            else
                cells(jj).area_crop_est = [];
            end
        end
    end
	% save with same filename
    save([pathname filename{hh}],'cells','-v7.3')
end

%% Estimate number of nucleoids for each frame
% Estimates this based on HU labelled area using Otsu's thresholding.
% Limits of nucleoid area are extracted from population experiments.

min_nuc_area = 50;      % minimum area of nucleoid in pixels
max_nuc_area = 800;     % maximum area of a single nucleoid in pixels

% get cells file
[filename,pathname] = uigetfile('*.mat', 'Select','Select','MultiSelect', 'on');
if ~iscell(filename)
    name = filename;
    filename = cell(1);
    filename{1} = name;
end

for hh = 1:length(filename)   
    load([pathname filename{hh}])
	for ii = 1:length(cells) % cells
        if isempty(cells(ii).area_crop_est) || (~isempty(cells(ii).area_crop_est) && length(cells(ii).area_crop_est) > 1)
            % pre-allocate field
            cells(ii).no_chrom_frame = NaN(1,length(cells(ii).fluo1));
            for jj = 1:length(cells(ii).fluo1) % frames
                if ~isempty(cells(ii).fluo1{jj})
                    % segment nucleoid using otsu
                    cell2thresh = immultiply(cells(ii).fluo1{jj},cells(ii).mask{jj});
                    thresh = multithresh(cell2thresh,2);
                    mask = cell2thresh > max(thresh);
                    % remove too small areas
                    mask2 = bwareaopen(mask, min_nuc_area);
                    % unique labels
                    L = bwlabel(mask2);
                    % chrom areas
                    areas = zeros(1,(length(unique(L))-1));
                    for kk = 1:(length(unique(L))-1)
                        areas(kk) = sum(L(:) == kk);
                    end

                    % save info to cells if areas smaller than max area
                    if sum(areas > max_nuc_area) == 0
                        cells(ii).no_chrom_frame(jj) = length(unique(L))-1;
                    else
                        cells(ii).no_chrom_frame(jj) = Inf;
                    end
                end
            end
            no_chrom_frame = cells(ii).no_chrom_frame;
            
            % define cropped version (if no crop then identical)
            if ~isempty(cells(ii).area_crop_est)
                cropped_no_chrom = no_chrom_frame(1:length(cells(ii).area_crop_est)); 
            else
                cropped_no_chrom = no_chrom_frame;
            end
            % remove NaNs
            cropped_no_chrom = cropped_no_chrom(~isnan(cropped_no_chrom));
            no_chrom_frame = no_chrom_frame(~isnan(no_chrom_frame));
            
            % replace Inf (too large areas) by the consequent value
            % this value must be more than 1 be definition
            % loop from end to beginning
            for kk = (length(no_chrom_frame)-1):-1:1
                if no_chrom_frame(kk) == Inf
                    no_chrom_frame(kk) = no_chrom_frame(kk+1);
                    % do the same for cropped if long enough
                    if kk <= length(cropped_no_chrom)
                        cropped_no_chrom(kk) = no_chrom_frame(kk+1);
                    end
                end 
            end
            
            % determine what is the maximum number of chromosomes from
            % cropped vector (Inf values are not valid)
            if isempty(cropped_no_chrom)
                cells(ii).no_chrom_est = NaN;
            elseif sum(cropped_no_chrom == Inf) > 0
                cells(ii).no_chrom_est = NaN;
            else
                cells(ii).no_chrom_est = cropped_no_chrom(end);
            end
  
        end
	end
    % save with same filename
    save([pathname filename{hh}],'cells','-v7.3')
end

%% Estimating growth rate and plotting Data for 1N cells

no_repeats = 3;             % number of repeats
time_int = 5;               % time interval (min)
win_size = 5;               % smoothing window size

limits = [2.2 10];          % plotting limits
interval = 0.5;             % interval for binning data
range = (limits(1)-interval/2):interval:(limits(2)+interval*1.5);
area_bins = range(1:(end-1))+interval/2;

step_mean = NaN(length(range)-1, no_repeats);
step_std = NaN(length(range)-1, no_repeats);
step_sem = NaN(length(range)-1, no_repeats);
rel_step_mean = NaN(length(range)-1, no_repeats);
rel_step_std = NaN(length(range)-1, no_repeats);
rel_step_sem = NaN(length(range)-1, no_repeats);
cell_size = cell(no_repeats,10000);
raw_size = cell(no_repeats,10000);
step_size = cell(no_repeats,10000);
rel_step_size = cell(no_repeats,10000);
IDs = cell(no_repeats,10000);
no_spots = cell(no_repeats,10000);
no_cells = [];
for ww = 1:no_repeats
    % get cells file
    [filename,pathname] = uigetfile('*.mat', 'Select','Select','MultiSelect', 'on');
    if ~iscell(filename)
        name = filename;
        filename = cell(1);
        filename{1} = name;
    end

    no_cells_single = 0;
    counter = 1;
    for hh = 1:length(filename)   
        load([pathname filename{hh}])

        for ii = 1:length(cells) % cells
            % find number of chromosomes at the max range (cropped area)
            if isfield(cells(ii), 'area_crop_est') && ~isempty(cells(ii).area_crop_est)
                areas = cells(ii).area_crop_est;
            else
                areas = cells(ii).area;
            end
            if isfield(cells(ii), 'no_chrom_frame')
                if max(areas) > max(range) && min(areas) < max(range)
                    inds = areas > max(range);
                    chrom_numbers = cells(ii).no_chrom_frame(inds);
                    chrom_numbers = chrom_numbers(~isnan(chrom_numbers));
                    if ~isempty(chrom_numbers)
                        no_chrom_range = chrom_numbers(1);
                    else
                        no_chrom_range = NaN;
                    end
                elseif min(areas) > max(range) % cell always too large for max range
                    no_chrom_range = NaN;
                else
                    no_chrom_range = cells(ii).no_chrom_est;
                end 
            end
                
            if ~isempty(cells(ii).no_chrom_est) && cells(ii).no_chrom_est == 1
                % spots field, replace NaN values by previous number
                if isfield(cells(ii), 'no_spots')
                    for kk = 1:length(cells(ii).no_spots)
                        if kk == 1 & isnan(cells(ii).no_spots(kk))
                            values = cells(ii).no_spots(~isnan(cells(ii).no_spots));
                            cells(ii).no_spots(kk) = values(1);
                        elseif isnan(cells(ii).no_spots(kk))
                            cells(ii).no_spots(kk) = cells(ii).no_spots(kk-1);
                        end
                    end
                end
                % cropped cells
                if isfield(cells,'area_crop_est') && ~isempty(cells(ii).area_crop_est)
                    sizes = cells(ii).area_crop_est;
                else
                    sizes = cells(ii).area;
                end
                
                % smooth area curve
                smoothed_sizes = smooth(sizes, win_size)';
                % truncate the ends by the window size
                smoothed_sizes = smoothed_sizes((floor(win_size/2)+1):(end-floor(win_size/2)));
                diff_size = diff(smoothed_sizes)./time_int;
                
                raw_size{ww,counter} = sizes;
                cell_size{ww,counter} = smoothed_sizes(2:end);
                step_size{ww,counter} = diff_size;
                rel_step_size{ww,counter} = (diff_size./smoothed_sizes(2:end));
                IDs{ww,counter} = ii;

                if isfield(cells(ii), 'no_spots')
                    values = cells(ii).no_spots(1:length(sizes));
                    % truncate the ends by the window size
                    values = values((floor(win_size/2)+1):(end-floor(win_size/2)));
                    no_spots{ww,counter} = values(1:(end-1));
                end
                no_cells_single = no_cells_single + 1;
                counter = counter + 1;
            end
        end
    end
    % which variable data is binned by
    xval = [cell_size{ww,:}];
    step_size_single = [step_size{ww,:}];
    rel_step_size_single = [rel_step_size{ww,:}];

    % bin step_size data and calculate mean, std, sem
    [B,~,idx] = histcounts(xval,range);
    means = accumarray(idx(idx > 0)',step_size_single(idx > 0),[],@mean);
    stds = accumarray(idx(idx > 0)',step_size_single(idx > 0),[],@std);
    sems = accumarray(idx(idx > 0)',step_size_single(idx > 0),[],@(x) std(x)./sqrt(length(x)));
	step_mean(B>0,ww) = means(means~=0);
    step_std(B>1,ww) = stds(stds~=0);
    step_sem(B>1,ww) = sems(sems~=0);
    
    % bin rel_step_size data and calculate mean, std, sem
    [B,~,idx] = histcounts(xval,range);
    means = accumarray(idx(idx > 0)',rel_step_size_single(idx > 0),[],@mean);
    stds = accumarray(idx(idx > 0)',rel_step_size_single(idx > 0),[],@std);
    sems = accumarray(idx(idx > 0)',rel_step_size_single(idx > 0),[],@(x) std(x)./sqrt(length(x)));
	rel_step_mean(B>0,ww) = means(means~=0);
    rel_step_std(B>1,ww) = stds(stds~=0);
    rel_step_sem(B>1,ww) = sems(sems~=0);

    % combine cell numbers
    no_cells = [no_cells no_cells_single];
end
% calculate mean, std and sem of means
step_size_mean_all = mean(step_mean,2);
step_size_std_all = std(step_mean,[],2);
step_size_sem_all = std(step_mean,[],2)./sqrt(no_repeats);
rel_step_size_mean_all = mean(rel_step_mean,2);
rel_step_size_std_all = std(rel_step_mean,[],2);
rel_step_size_sem_all = std(rel_step_mean,[],2)./sqrt(no_repeats);

% plotting
% absolute growth - SD
figure;
hold on;plot(area_bins,step_size_mean_all,'k','LineWidth',1)
hold on;plot(area_bins,step_size_mean_all-step_size_std_all,'k','LineWidth',0.5)
hold on;plot(area_bins,step_size_mean_all+step_size_std_all,'k','LineWidth',0.5)
% absolute growth - SEM
hold on;plot(area_bins,step_size_mean_all,'k','LineWidth',1)
hold on;plot(area_bins,step_size_mean_all-step_size_sem_all,'k','LineWidth',0.5)
hold on;plot(area_bins,step_size_mean_all+step_size_sem_all,'k','LineWidth',0.5)
% relative growth - SD
hold on;plot(area_bins,rel_step_size_mean_all,'k','LineWidth',1)
hold on;plot(area_bins,rel_step_size_mean_all-rel_step_size_std_all,'k','LineWidth',0.5)
hold on;plot(area_bins,rel_step_size_mean_all+rel_step_size_std_all,'k','LineWidth',0.5)
% relative growth - SEM
hold on;plot(area_bins,rel_step_size_mean_all,'k','LineWidth',1)
hold on;plot(area_bins,rel_step_size_mean_all-rel_step_size_sem_all,'k','LineWidth',0.5)
hold on;plot(area_bins,rel_step_size_mean_all+rel_step_size_sem_all,'k','LineWidth',0.5)

% add line for exponential growth based on normalized growth rate (exp growth)
growth_rate = mean([rel_step_size{:}]);
% calculate absolute growth and plot
y = growth_rate.*area_bins;
figure;plot(area_bins,y,'r')
figure;plot([area_bins(1) area_bins(end)],[growth_rate growth_rate],'r')

%check quantiles
quantiles = quantile([cell_size{:}],[0.025 0.975]);

xval = [cell_size{:}];
% xval = times;

yval = [step_size{:}];
% yval = [rel_step_size{:}];

% plot scatter
figure;plot(xval,yval,'.','Color',[0.8 0.8 0.8],'MarkerSize',12)
xlabel('Cell area (\mum^2)')
ylabel('Growth rate (\mum^2/min)')
title([num2str(sum(no_cells)) ' cells'])
axis([1 10 0 Inf])

% ksdensity for scatter
figure;
density = ksdensity([xval',yval'], [xval',yval']);
scatter(xval', yval', 20, density,'filled');
xlabel('Cell area (\mum^2)')
ylabel('Growth rate (\mum^2/min)')
title([num2str(sum(no_cells)) ' cells'])
axis([1 10 0 Inf])

% plot growth for each cell
figure; hold on;
for ii = 1:size(cell_size,1)
    for jj = 1:size(cell_size,2)
        if ~isempty(cell_size{ii,jj})
            plot(cell_size{ii,jj},step_size{ii,jj},'Color',[0.8 0.8 0.8],'LineWidth',1,'DisplayName',[num2str(ii) '-' num2str(IDs{ii,jj})])
        end
    end
end
title([num2str(sum(no_cells)) ' cells'])
axis([1 10 0 Inf])

% relative linear growth
hold on; plot(area_bins,area_bins(1)*0.0161628./area_bins ,'k','LineWidth',1)

% plot cell areas
figure; hold on;
for ii = 1:size(cell_size,1)
    for jj = 1:size(cell_size,2)
        if ~isempty(cell_size{ii,jj})
            plot(time_int.*(0:(length(raw_size{ii,jj}))-1),raw_size{ii,jj},'b','LineWidth',1,'DisplayName',[num2str(ii) '-' num2str(IDs{ii,jj})])
        end
    end
end
title([num2str(sum(no_cells)) ' cells'])

%% Estimating growth rate and plotting Data for WT cells

no_repeats = 3;             % number of repeats
time_int = 5;               % time interval (min)
win_size = 5;               % smoothing window size

limits = [2.2 10];          % plotting limits
interval = 0.5;             % interval for binning data
range = (limits(1)-interval/2):interval:(limits(2)+interval*1.5);
area_bins = range(1:(end-1))+interval/2;

step_mean = NaN(length(range)-1, no_repeats);
step_std = NaN(length(range)-1, no_repeats);
step_sem = NaN(length(range)-1, no_repeats);
rel_step_mean = NaN(length(range)-1, no_repeats);
rel_step_std = NaN(length(range)-1, no_repeats);
rel_step_sem = NaN(length(range)-1, no_repeats);
cell_size = cell(no_repeats,10000);
raw_size = cell(no_repeats,10000);
step_size = cell(no_repeats,10000);
rel_step_size = cell(no_repeats,10000);
IDs = cell(no_repeats,10000);
no_spots = cell(no_repeats,10000);
no_cells = [];
for ww = 1:no_repeats
    % get cells file
    [filename,pathname] = uigetfile('*.mat', 'Select','Select','MultiSelect', 'on');
    if ~iscell(filename)
        name = filename;
        filename = cell(1);
        filename{1} = name;
    end

    no_cells_single = 0;
    counter = 1;
    for hh = 1:length(filename)   
        load([pathname filename{hh}])

        for ii = 1:length(cells) % cells
            % find number of chromosomes at the max range (cropped area)
            if isfield(cells(ii), 'area_crop_est') && ~isempty(cells(ii).area_crop_est)
                areas = cells(ii).area_crop_est;
            else
                areas = cells(ii).area;
            end
            if isfield(cells(ii), 'no_chrom_frame')
                if max(areas) > max(range) && min(areas) < max(range)
                    inds = areas > max(range);
                    chrom_numbers = cells(ii).no_chrom_frame(inds);
                    chrom_numbers = chrom_numbers(~isnan(chrom_numbers));
                    if ~isempty(chrom_numbers)
                        no_chrom_range = chrom_numbers(1);
                    else
                        no_chrom_range = NaN;
                    end
                elseif min(areas) > max(range) % cell always too large for max range
                    no_chrom_range = NaN;
                else
                    no_chrom_range = cells(ii).no_chrom_est;
                end 
            end
                
            if length(cells(ii).daughters) == 2 || (~isempty(cells(ii).parent) && cells(ii).birth > 1)
                % spots field, replace NaN values by previous number
                if isfield(cells(ii), 'no_spots')
                    for kk = 1:length(cells(ii).no_spots)
                        if kk == 1 & isnan(cells(ii).no_spots(kk))
                            values = cells(ii).no_spots(~isnan(cells(ii).no_spots));
                            cells(ii).no_spots(kk) = values(1);
                        elseif isnan(cells(ii).no_spots(kk))
                            cells(ii).no_spots(kk) = cells(ii).no_spots(kk-1);
                        end
                    end
                end
                % cropped cells
                if isfield(cells,'area_crop_est') && ~isempty(cells(ii).area_crop_est)
                    sizes = cells(ii).area_crop_est;
                else
                    sizes = cells(ii).area;
                end
                
                % smooth area curve
                smoothed_sizes = smooth(sizes, win_size)';
                % truncate the ends by the window size
                smoothed_sizes = smoothed_sizes((floor(win_size/2)+1):(end-floor(win_size/2)));
                diff_size = diff(smoothed_sizes)./time_int;
                
                raw_size{ww,counter} = sizes;
                cell_size{ww,counter} = smoothed_sizes(2:end);
                step_size{ww,counter} = diff_size;
                rel_step_size{ww,counter} = (diff_size./smoothed_sizes(2:end));
                IDs{ww,counter} = ii;

                if isfield(cells(ii), 'no_spots')
                    values = cells(ii).no_spots(1:length(sizes));
                    % truncate the ends by the window size
                    values = values((floor(win_size/2)+1):(end-floor(win_size/2)));
                    no_spots{ww,counter} = values(1:(end-1));
                end
                no_cells_single = no_cells_single + 1;
                counter = counter + 1;
            end
        end
    end
    % which variable data is binned by
    xval = [cell_size{ww,:}];
    step_size_single = [step_size{ww,:}];
    rel_step_size_single = [rel_step_size{ww,:}];

    % bin step_size data and calculate mean, std, sem
    [B,~,idx] = histcounts(xval,range);
    means = accumarray(idx(idx > 0)',step_size_single(idx > 0),[],@mean);
    stds = accumarray(idx(idx > 0)',step_size_single(idx > 0),[],@std);
    sems = accumarray(idx(idx > 0)',step_size_single(idx > 0),[],@(x) std(x)./sqrt(length(x)));
	step_mean(B>0,ww) = means(means~=0);
    step_std(B>1,ww) = stds(stds~=0);
    step_sem(B>1,ww) = sems(sems~=0);
    
    % bin rel_step_size data and calculate mean, std, sem
    [B,~,idx] = histcounts(xval,range);
    means = accumarray(idx(idx > 0)',rel_step_size_single(idx > 0),[],@mean);
    stds = accumarray(idx(idx > 0)',rel_step_size_single(idx > 0),[],@std);
    sems = accumarray(idx(idx > 0)',rel_step_size_single(idx > 0),[],@(x) std(x)./sqrt(length(x)));
	rel_step_mean(B>0,ww) = means(means~=0);
    rel_step_std(B>1,ww) = stds(stds~=0);
    rel_step_sem(B>1,ww) = sems(sems~=0);

    % combine cell numbers
    no_cells = [no_cells no_cells_single];
end
% calculate mean, std and sem of means
step_size_mean_all = mean(step_mean,2);
step_size_std_all = std(step_mean,[],2);
step_size_sem_all = std(step_mean,[],2)./sqrt(no_repeats);
rel_step_size_mean_all = mean(rel_step_mean,2);
rel_step_size_std_all = std(rel_step_mean,[],2);
rel_step_size_sem_all = std(rel_step_mean,[],2)./sqrt(no_repeats);

% plotting
% absolute growth - SD
figure;
hold on;plot(area_bins,step_size_mean_all,'k','LineWidth',1)
hold on;plot(area_bins,step_size_mean_all-step_size_std_all,'k','LineWidth',0.5)
hold on;plot(area_bins,step_size_mean_all+step_size_std_all,'k','LineWidth',0.5)
% absolute growth - SEM
hold on;plot(area_bins,step_size_mean_all,'k','LineWidth',1)
hold on;plot(area_bins,step_size_mean_all-step_size_sem_all,'k','LineWidth',0.5)
hold on;plot(area_bins,step_size_mean_all+step_size_sem_all,'k','LineWidth',0.5)
% relative growth - SD
hold on;plot(area_bins,rel_step_size_mean_all,'k','LineWidth',1)
hold on;plot(area_bins,rel_step_size_mean_all-rel_step_size_std_all,'k','LineWidth',0.5)
hold on;plot(area_bins,rel_step_size_mean_all+rel_step_size_std_all,'k','LineWidth',0.5)
% relative growth - SEM
hold on;plot(area_bins,rel_step_size_mean_all,'k','LineWidth',1)
hold on;plot(area_bins,rel_step_size_mean_all-rel_step_size_sem_all,'k','LineWidth',0.5)
hold on;plot(area_bins,rel_step_size_mean_all+rel_step_size_sem_all,'k','LineWidth',0.5)

% add line for exponential growth based on normalized growth rate (exp growth)
growth_rate = mean([rel_step_size{:}]);
% calculate absolute growth and plot
y = growth_rate.*area_bins;
figure;plot(area_bins,y,'r')
figure;plot([area_bins(1) area_bins(end)],[growth_rate growth_rate],'r')

%check quantiles
quantiles = quantile([cell_size{:}],[0.025 0.975]);

xval = [cell_size{:}];
% xval = times;

yval = [step_size{:}];
% yval = [rel_step_size{:}];

% plot scatter
figure;plot(xval,yval,'.','Color',[0.8 0.8 0.8],'MarkerSize',12)
xlabel('Cell area (\mum^2)')
ylabel('Growth rate (\mum^2/min)')
title([num2str(sum(no_cells)) ' cells'])
axis([1 10 0 Inf])

% ksdensity for scatter
figure;
density = ksdensity([xval',yval'], [xval',yval']);
scatter(xval', yval', 20, density,'filled');
xlabel('Cell area (\mum^2)')
ylabel('Growth rate (\mum^2/min)')
title([num2str(sum(no_cells)) ' cells'])
axis([1 10 0 Inf])

% plot growth for each cell
figure; hold on;
for ii = 1:size(cell_size,1)
    for jj = 1:size(cell_size,2)
        if ~isempty(cell_size{ii,jj})
            plot(cell_size{ii,jj},step_size{ii,jj},'Color',[0.8 0.8 0.8],'LineWidth',1,'DisplayName',[num2str(ii) '-' num2str(IDs{ii,jj})])
        end
    end
end
title([num2str(sum(no_cells)) ' cells'])
axis([1 10 0 Inf])

% relative linear growth
hold on; plot(area_bins,area_bins(1)*0.0161628./area_bins ,'k','LineWidth',1)

% plot cell areas
figure; hold on;
for ii = 1:size(cell_size,1)
    for jj = 1:size(cell_size,2)
        if ~isempty(cell_size{ii,jj})
            plot(time_int.*(0:(length(raw_size{ii,jj}))-1),raw_size{ii,jj},'b','LineWidth',1,'DisplayName',[num2str(ii) '-' num2str(IDs{ii,jj})])
        end
    end
end
title([num2str(sum(no_cells)) ' cells'])

%% Estimating growth rate and plotting data for dnaC2 cells

no_repeats = 3;             % number of repeats
time_int = 5;               % time interval (min)
win_size = 5;               % smoothing window size

limits = [2.2 10];          % plotting limits
interval = 0.5;             % interval for binning data
range = (limits(1)-interval/2):interval:(limits(2)+interval*1.5);
area_bins = range(1:(end-1))+interval/2;

step_mean = NaN(length(range)-1, no_repeats);
step_std = NaN(length(range)-1, no_repeats);
step_sem = NaN(length(range)-1, no_repeats);
rel_step_mean = NaN(length(range)-1, no_repeats);
rel_step_std = NaN(length(range)-1, no_repeats);
rel_step_sem = NaN(length(range)-1, no_repeats);
cell_size = cell(no_repeats,10000);
raw_size = cell(no_repeats,10000);
step_size = cell(no_repeats,10000);
rel_step_size = cell(no_repeats,10000);
IDs = cell(no_repeats,10000);
no_spots = cell(no_repeats,10000);
no_cells = [];
for ww = 1:no_repeats
    % get cells file
    [filename,pathname] = uigetfile('*.mat', 'Select','Select','MultiSelect', 'on');
    if ~iscell(filename)
        name = filename;
        filename = cell(1);
        filename{1} = name;
    end

    no_cells_single = 0;
    counter = 1;
    for hh = 1:length(filename)   
        load([pathname filename{hh}])

        for ii = 1:length(cells) % cells
            % find number of chromosomes at the max range (cropped area)
            if isfield(cells(ii), 'area_crop_est') && ~isempty(cells(ii).area_crop_est)
                areas = cells(ii).area_crop_est;
            else
                areas = cells(ii).area;
            end
            if isfield(cells(ii), 'no_chrom_frame')
                if max(areas) > max(range) && min(areas) < max(range)
                    inds = areas > max(range);
                    chrom_numbers = cells(ii).no_chrom_frame(inds);
                    chrom_numbers = chrom_numbers(~isnan(chrom_numbers));
                    if ~isempty(chrom_numbers)
                        no_chrom_range = chrom_numbers(1);
                    else
                        no_chrom_range = NaN;
                    end
                elseif min(areas) > max(range) % cell always too large for max range
                    no_chrom_range = NaN;
                else
                    no_chrom_range = cells(ii).no_chrom_est;
                end 
            end

            if ~isnan(no_chrom_range) && ~isinf(no_chrom_range) && no_chrom_range > 2
                % spots field, replace NaN values by previous number
                if isfield(cells(ii), 'no_spots')
                    for kk = 1:length(cells(ii).no_spots)
                        if kk == 1 & isnan(cells(ii).no_spots(kk))
                            values = cells(ii).no_spots(~isnan(cells(ii).no_spots));
                            cells(ii).no_spots(kk) = values(1);
                        elseif isnan(cells(ii).no_spots(kk))
                            cells(ii).no_spots(kk) = cells(ii).no_spots(kk-1);
                        end
                    end
                end
                % cropped cells
                if isfield(cells,'area_crop_est') && ~isempty(cells(ii).area_crop_est)
                    sizes = cells(ii).area_crop_est;
                else
                    sizes = cells(ii).area;
                end
                
                % smooth area curve
                smoothed_sizes = smooth(sizes, win_size)';
                % truncate the ends by the window size
                smoothed_sizes = smoothed_sizes((floor(win_size/2)+1):(end-floor(win_size/2)));
                diff_size = diff(smoothed_sizes)./time_int;
                
                raw_size{ww,counter} = sizes;
                cell_size{ww,counter} = smoothed_sizes(2:end);
                step_size{ww,counter} = diff_size;
                rel_step_size{ww,counter} = (diff_size./smoothed_sizes(2:end));
                IDs{ww,counter} = ii;

                if isfield(cells(ii), 'no_spots')
                    values = cells(ii).no_spots(1:length(sizes));
                    % truncate the ends by the window size
                    values = values((floor(win_size/2)+1):(end-floor(win_size/2)));
                    no_spots{ww,counter} = values(1:(end-1));
                end
                no_cells_single = no_cells_single + 1;
                counter = counter + 1;
            end
        end
    end
    % which variable data is binned by
    xval = [cell_size{ww,:}];
    step_size_single = [step_size{ww,:}];
    rel_step_size_single = [rel_step_size{ww,:}];

    % bin step_size data and calculate mean, std, sem
    [B,~,idx] = histcounts(xval,range);
    means = accumarray(idx(idx > 0)',step_size_single(idx > 0),[],@mean);
    stds = accumarray(idx(idx > 0)',step_size_single(idx > 0),[],@std);
    sems = accumarray(idx(idx > 0)',step_size_single(idx > 0),[],@(x) std(x)./sqrt(length(x)));
	step_mean(B>0,ww) = means(means~=0);
    step_std(B>1,ww) = stds(stds~=0);
    step_sem(B>1,ww) = sems(sems~=0);
    
    % bin rel_step_size data and calculate mean, std, sem
    [B,~,idx] = histcounts(xval,range);
    means = accumarray(idx(idx > 0)',rel_step_size_single(idx > 0),[],@mean);
    stds = accumarray(idx(idx > 0)',rel_step_size_single(idx > 0),[],@std);
    sems = accumarray(idx(idx > 0)',rel_step_size_single(idx > 0),[],@(x) std(x)./sqrt(length(x)));
	rel_step_mean(B>0,ww) = means(means~=0);
    rel_step_std(B>1,ww) = stds(stds~=0);
    rel_step_sem(B>1,ww) = sems(sems~=0);

    % combine cell numbers
    no_cells = [no_cells no_cells_single];
end
% calculate mean, std and sem of means
step_size_mean_all = mean(step_mean,2);
step_size_std_all = std(step_mean,[],2);
step_size_sem_all = std(step_mean,[],2)./sqrt(no_repeats);
rel_step_size_mean_all = mean(rel_step_mean,2);
rel_step_size_std_all = std(rel_step_mean,[],2);
rel_step_size_sem_all = std(rel_step_mean,[],2)./sqrt(no_repeats);

% plotting
% absolute growth - SD
figure;
hold on;plot(area_bins,step_size_mean_all,'k','LineWidth',1)
hold on;plot(area_bins,step_size_mean_all-step_size_std_all,'k','LineWidth',0.5)
hold on;plot(area_bins,step_size_mean_all+step_size_std_all,'k','LineWidth',0.5)
% absolute growth - SEM
hold on;plot(area_bins,step_size_mean_all,'k','LineWidth',1)
hold on;plot(area_bins,step_size_mean_all-step_size_sem_all,'k','LineWidth',0.5)
hold on;plot(area_bins,step_size_mean_all+step_size_sem_all,'k','LineWidth',0.5)
% relative growth - SD
hold on;plot(area_bins,rel_step_size_mean_all,'k','LineWidth',1)
hold on;plot(area_bins,rel_step_size_mean_all-rel_step_size_std_all,'k','LineWidth',0.5)
hold on;plot(area_bins,rel_step_size_mean_all+rel_step_size_std_all,'k','LineWidth',0.5)
% relative growth - SEM
hold on;plot(area_bins,rel_step_size_mean_all,'k','LineWidth',1)
hold on;plot(area_bins,rel_step_size_mean_all-rel_step_size_sem_all,'k','LineWidth',0.5)
hold on;plot(area_bins,rel_step_size_mean_all+rel_step_size_sem_all,'k','LineWidth',0.5)

% add line for exponential growth based on normalized growth rate (exp growth)
growth_rate = mean([rel_step_size{:}]);
% calculate absolute growth and plot
y = growth_rate.*area_bins;
figure;plot(area_bins,y,'r')
figure;plot([area_bins(1) area_bins(end)],[growth_rate growth_rate],'r')

%check quantiles
quantiles = quantile([cell_size{:}],[0.025 0.975]);

xval = [cell_size{:}];
% xval = times;

yval = [step_size{:}];
% yval = [rel_step_size{:}];

% plot scatter
figure;plot(xval,yval,'.','Color',[0.8 0.8 0.8],'MarkerSize',12)
xlabel('Cell area (\mum^2)')
ylabel('Growth rate (\mum^2/min)')
title([num2str(sum(no_cells)) ' cells'])
axis([1 10 0 Inf])

% ksdensity for scatter
figure;
density = ksdensity([xval',yval'], [xval',yval']);
scatter(xval', yval', 20, density,'filled');
xlabel('Cell area (\mum^2)')
ylabel('Growth rate (\mum^2/min)')
title([num2str(sum(no_cells)) ' cells'])
axis([1 10 0 Inf])

% plot growth for each cell
figure; hold on;
for ii = 1:size(cell_size,1)
    for jj = 1:size(cell_size,2)
        if ~isempty(cell_size{ii,jj})
            plot(cell_size{ii,jj},step_size{ii,jj},'Color',[0.8 0.8 0.8],'LineWidth',1,'DisplayName',[num2str(ii) '-' num2str(IDs{ii,jj})])
        end
    end
end
title([num2str(sum(no_cells)) ' cells'])
axis([1 10 0 Inf])

% relative linear growth
hold on; plot(area_bins,area_bins(1)*0.0161628./area_bins ,'k','LineWidth',1)

% plot cell areas
figure; hold on;
for ii = 1:size(cell_size,1)
    for jj = 1:size(cell_size,2)
        if ~isempty(cell_size{ii,jj})
            plot(time_int.*(0:(length(raw_size{ii,jj}))-1),raw_size{ii,jj},'b','LineWidth',1,'DisplayName',[num2str(ii) '-' num2str(IDs{ii,jj})])
        end
    end
end
title([num2str(sum(no_cells)) ' cells'])

%% Check for constant growth in asynchronously dividing cells
% calculates normalized growth rate at different frames

% get file
[filename,pathname] = uigetfile('*info.mat', 'Select','Select');
load([pathname filename])

norm_rates = [];
frames = [];
for ii = 1:length(cell_info)
    if length(cell_info(ii).area) > 1
        rate = diff(cell_info(ii).area);
        norm_rate = rate./cell_info(ii).area(1:(end-1));
        frame = cell_info(ii).birth:(cell_info(ii).division-1);
        norm_rates = [norm_rates norm_rate];
        frames = [frames frame];
    end
end

% bin the data
[B,idx] = histc(frames,1:max([cell_info.division]));
Y_mean = accumarray(idx(:),norm_rates,[],@median);
X_mean = accumarray(idx(:),frames,[],@median);

figure;hold on;
plot(frames,norm_rates,'.k')
plot(X_mean,Y_mean,'r','LineWidth',2)
axis([1 max(X_mean) 0 2*max(Y_mean)])
