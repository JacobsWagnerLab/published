%% Time-lapse analysis using cell volumes
% Use this script to convert cell areas to volume using Oufti script.
% Run the normal area based analysis before this (MAIN_timelapse.m).

%% Estimate volume from cell mask using Oufti
% this script adds a new field into cells structure 'volume_est'

pixel = 0.0666;     % pixel size in um

% Oufti parameters required
STEP_SIZE = 1.2;
TOLERANCE = 1e-5;
MESH_WIDTH = 500;

% select *cells.mat data file
[filenameC,pathnameC] = uigetfile('*.mat', 'Select cells file','MultiSelect', 'off');
load([pathnameC filenameC])

disp(['Total ' num2str(length(cells))])
% loop through individual cells
if ~isfield(cells, 'volume_est')
    cells(1).volume_est = [];
end
for ii = 1:length(cells)
    % only cells that don't have volume already
    if isempty(cells(ii).volume_est)
        disp(num2str(ii))

        % loop through time points
        if ~isempty(cells(ii).area_crop_est)
            max_frame = length(cells(ii).area_crop_est);
        else
            max_frame = length(cells(ii).area);
        end
        for jj = 1:max_frame
            mask = cells(ii).mask{jj};

            % convert mask into Oufti mesh
            [meshes,models,unique_bits] = mask2mesh(mask,STEP_SIZE,TOLERANCE,MESH_WIDTH);
            if isempty(meshes)
                break
            else
                mesh = meshes{1,1};
                
                % estimate cell length
                steplength = edist(mesh(2:end,1)+mesh(2:end,3),mesh(2:end,2)+mesh(2:end,4),...
                    mesh(1:end-1,1)+mesh(1:end-1,3),mesh(1:end-1,2)+mesh(1:end-1,4))/2;
                tot_length = sum(steplength);

                % estimate cell volume
                d = edist(mesh(:,1),mesh(:,2),mesh(:,3),mesh(:,4));
                stepvolume = (d(1:end-1).*d(2:end) + (d(1:end-1)-d(2:end)).^2/3).*steplength*pi/4;
                volume_pixel = sum(stepvolume);

                % volume in um3
                volume = volume_pixel.*pixel^3;

                % save volume into cells structure
                cells(ii).volume_est(jj) = volume;
            end
        end
    end
end

% save file
save([pathnameC filenameC],'cells','-v7.3')

%% Estimating growth rate and plotting Data for 1N cells

no_repeats = 3;             % number of repeats
time_int = 5;               % time interval (min)
win_size = 5;               % smoothing window size

% binning
limits = [0.54 6.5];          % plotting limits
interval = 0.5;             % interval for binning data
range = (limits(1)-interval/2):interval:(limits(2)+interval*1.5);
volume_bins = range(1:(end-1))+interval/2;

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
            volumes = cells(ii).volume_est;
            
            if ~isempty(volumes)
                if isfield(cells(ii), 'no_chrom_frame')
                    if max(volumes) > max(range) && min(volumes) < max(range)
                        inds = volumes > max(range);
                        chrom_numbers = cells(ii).no_chrom_frame(inds);
                        chrom_numbers = chrom_numbers(~isnan(chrom_numbers));
                        if ~isempty(chrom_numbers)
                            no_chrom_range = chrom_numbers(1);
                        else
                            no_chrom_range = NaN;
                        end
                    elseif min(volumes) > max(range) % cell always too large for max range
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
                    sizes = cells(ii).volume_est;

                    % smooth curve
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
hold on;plot(volume_bins,step_size_mean_all,'r','LineWidth',1)
hold on;plot(volume_bins,step_size_mean_all-step_size_std_all,'r','LineWidth',0.5)
hold on;plot(volume_bins,step_size_mean_all+step_size_std_all,'r','LineWidth',0.5)
% absolute growth - SEM
% hold on;plot(volume_bins,step_size_mean_all,'k','LineWidth',1)
% hold on;plot(volume_bins,step_size_mean_all-step_size_sem_all,'k','LineWidth',0.5)
% hold on;plot(volume_bins,step_size_mean_all+step_size_sem_all,'k','LineWidth',0.5)
% relative growth - SD
hold on;plot(volume_bins,rel_step_size_mean_all,'r','LineWidth',1)
hold on;plot(volume_bins,rel_step_size_mean_all-rel_step_size_std_all,'r','LineWidth',0.5)
hold on;plot(volume_bins,rel_step_size_mean_all+rel_step_size_std_all,'r','LineWidth',0.5)
% relative growth - SEM
% hold on;plot(volume_bins,rel_step_size_mean_all,'k','LineWidth',1)
% hold on;plot(volume_bins,rel_step_size_mean_all-rel_step_size_sem_all,'k','LineWidth',0.5)
% hold on;plot(volume_bins,rel_step_size_mean_all+rel_step_size_sem_all,'k','LineWidth',0.5)

% add line for exponential growth based on normalized growth rate (exp growth)
growth_rate = mean([rel_step_size{:}]);
% calculate absolute growth and plot
y = growth_rate.*volume_bins;
figure;plot(volume_bins,y,'r')
figure;plot([volume_bins(1) volume_bins(end)],[growth_rate growth_rate],'r')

%check quantiles
quantiles = quantile([cell_size{:}],[0.025 0.975]);

xval = [cell_size{:}];
% xval = times;

yval = [step_size{:}];
% yval = [rel_step_size{:}];

% plot scatter
figure;plot(xval,yval,'.','Color',[0.8 0.8 0.8],'MarkerSize',12)
xlabel('Cell volume (\mum^3)')
ylabel('Growth rate (\mum^3/min)')
title([num2str(sum(no_cells)) ' cells'])
axis([1 10 0 Inf])

% ksdensity for scatter
figure;
density = ksdensity([xval',yval'], [xval',yval']);
scatter(xval', yval', 20, density,'filled');
xlabel('Cell volume (\mum^3)')
ylabel('Growth rate (\mum^3/min)')
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
hold on; plot(volume_bins,volume_bins(1)*0.0161628./volume_bins ,'k','LineWidth',1)

% plot cell volumes
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

limits = [0.54 1.83];          % plotting limits
interval = 0.2;             % interval for binning data
range = (limits(1)-interval/2):interval:(limits(2)+interval*1.5);
volume_bins = range(1:(end-1))+interval/2;

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
            volumes = cells(ii).volume_est;
                
            if isfield(cells(ii), 'no_chrom_frame')
                if max(volumes) > max(range) && min(volumes) < max(range)
                    inds = volumes > max(range);
                    chrom_numbers = cells(ii).no_chrom_frame(inds);
                    chrom_numbers = chrom_numbers(~isnan(chrom_numbers));
                    if ~isempty(chrom_numbers)
                        no_chrom_range = chrom_numbers(1);
                    else
                        no_chrom_range = NaN;
                    end
                elseif min(volumes) > max(range) % cell always too large for max range
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
                sizes = cells(ii).volume_est;
                
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
hold on;plot(volume_bins,step_size_mean_all,'k','LineWidth',1)
hold on;plot(volume_bins,step_size_mean_all-step_size_std_all,'k','LineWidth',0.5)
hold on;plot(volume_bins,step_size_mean_all+step_size_std_all,'k','LineWidth',0.5)
% absolute growth - SEM
% hold on;plot(volume_bins,step_size_mean_all,'k','LineWidth',1)
% hold on;plot(volume_bins,step_size_mean_all-step_size_sem_all,'k','LineWidth',0.5)
% hold on;plot(volume_bins,step_size_mean_all+step_size_sem_all,'k','LineWidth',0.5)
% relative growth - SD
hold on;plot(volume_bins,rel_step_size_mean_all,'k','LineWidth',1)
hold on;plot(volume_bins,rel_step_size_mean_all-rel_step_size_std_all,'k','LineWidth',0.5)
hold on;plot(volume_bins,rel_step_size_mean_all+rel_step_size_std_all,'k','LineWidth',0.5)
% relative growth - SEM
% hold on;plot(volume_bins,rel_step_size_mean_all,'k','LineWidth',1)
% hold on;plot(volume_bins,rel_step_size_mean_all-rel_step_size_sem_all,'k','LineWidth',0.5)
% hold on;plot(volume_bins,rel_step_size_mean_all+rel_step_size_sem_all,'k','LineWidth',0.5)

% add line for exponential growth based on normalized growth rate (exp growth)
growth_rate = mean([rel_step_size{:}]);
% calculate absolute growth and plot
y = growth_rate.*volume_bins;
figure;plot(volume_bins,y,'r')
figure;plot([volume_bins(1) volume_bins(end)],[growth_rate growth_rate],'r')

%check quantiles
quantiles = quantile([cell_size{:}],[0.025 0.975]);

xval = [cell_size{:}];
% xval = times;

yval = [step_size{:}];
% yval = [rel_step_size{:}];

% plot scatter
figure;plot(xval,yval,'.','Color',[0.8 0.8 0.8],'MarkerSize',12)
xlabel('Cell volume (\mum^3)')
ylabel('Growth rate (\mum^3/min)')
title([num2str(sum(no_cells)) ' cells'])
axis([1 10 0 Inf])

% ksdensity for scatter
figure;
density = ksdensity([xval',yval'], [xval',yval']);
scatter(xval', yval', 20, density,'filled');
xlabel('Cell volume (\mum^3)')
ylabel('Growth rate (\mum^3/min)')
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
hold on; plot(volume_bins,volume_bins(1)*0.0161628./volume_bins ,'k','LineWidth',1)

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

%% Estimating growth rate and plotting Data for dnaC2 cells

no_repeats = 3;             % number of repeats
time_int = 5;               % time interval (min)
win_size = 5;               % smoothing window size

limits = [1.28 6.5];          % plotting limits
interval = 0.5;             % interval for binning data
range = (limits(1)-interval/2):interval:(limits(2)+interval*1.5);
volume_bins = range(1:(end-1))+interval/2;

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
            volumes = cells(ii).volume_est;

            if isfield(cells(ii), 'no_chrom_frame')
                if max(volumes) > max(range) && min(volumes) < max(range)
                    inds = volumes > max(range);
                    chrom_numbers = cells(ii).no_chrom_frame(inds);
                    chrom_numbers = chrom_numbers(~isnan(chrom_numbers));
                    if ~isempty(chrom_numbers)
                        no_chrom_range = chrom_numbers(1);
                    else
                        no_chrom_range = NaN;
                    end
                elseif min(volumes) > max(range) % cell always too large for max range
                    no_chrom_range = NaN;
                else
                    no_chrom_range = cells(ii).no_chrom_est;
                end 
            end

            %%%% change here the condition for no_chrom_range to match N
            if ~isnan(no_chrom_range) && ~isinf(no_chrom_range) && no_chrom_range == 1
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
                sizes = cells(ii).volume_est;
                
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
hold on;plot(volume_bins,step_size_mean_all,'k','LineWidth',1)
hold on;plot(volume_bins,step_size_mean_all-step_size_std_all,'k','LineWidth',0.5)
hold on;plot(volume_bins,step_size_mean_all+step_size_std_all,'k','LineWidth',0.5)
% absolute growth - SEM
hold on;plot(volume_bins,step_size_mean_all,'k','LineWidth',1)
hold on;plot(volume_bins,step_size_mean_all-step_size_sem_all,'k','LineWidth',0.5)
hold on;plot(volume_bins,step_size_mean_all+step_size_sem_all,'k','LineWidth',0.5)
% relative growth - SD
hold on;plot(volume_bins,rel_step_size_mean_all,'k','LineWidth',1)
hold on;plot(volume_bins,rel_step_size_mean_all-rel_step_size_std_all,'k','LineWidth',0.5)
hold on;plot(volume_bins,rel_step_size_mean_all+rel_step_size_std_all,'k','LineWidth',0.5)
% relative growth - SEM
hold on;plot(volume_bins,rel_step_size_mean_all,'k','LineWidth',1)
hold on;plot(volume_bins,rel_step_size_mean_all-rel_step_size_sem_all,'k','LineWidth',0.5)
hold on;plot(volume_bins,rel_step_size_mean_all+rel_step_size_sem_all,'k','LineWidth',0.5)

% add line for exponential growth based on normalized growth rate (exp growth)
growth_rate = mean([rel_step_size{:}]);
% calculate absolute growth and plot
y = growth_rate.*volume_bins;
figure;plot(volume_bins,y,'r')
figure;plot([volume_bins(1) volume_bins(end)],[growth_rate growth_rate],'r')

%check quantiles
quantiles = quantile([cell_size{:}],[0.025 0.975]);

xval = [cell_size{:}];
% xval = times;

yval = [step_size{:}];
% yval = [rel_step_size{:}];

% plot scatter
figure;plot(xval,yval,'.','Color',[0.8 0.8 0.8],'MarkerSize',12)
xlabel('Cell volume (\mum^3)')
ylabel('Growth rate (\mum^3/min)')
title([num2str(sum(no_cells)) ' cells'])
axis([1 10 0 Inf])

% ksdensity for scatter
figure;
density = ksdensity([xval',yval'], [xval',yval']);
scatter(xval', yval', 20, density,'filled');
xlabel('Cell volume (\mum^3)')
ylabel('Growth rate (\mum^3/min)')
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
hold on; plot(volume_bins,volume_bins(1)*0.0161628./volume_bins ,'k','LineWidth',1)

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
