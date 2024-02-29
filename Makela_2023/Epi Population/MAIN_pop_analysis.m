%% Population level analysis
% Detect cell areas using U-net (or some other algorithm)
% In each frame, each cell mask should have a unique mask ID.

% Assumes the following filename structure: 
%  ***xy**c*.tif
%       *c1.tif - phase contrast
%       *c1_mask.tif - mask file with unique mask numbers
%       *c2.tif - fluo1 channel
%       *c3.itf - fluo2 channel
%       ...

% 7/7/22 Jarno Makela

%% Collect cell based image data
% Save data into a .mat file. See options of the function

collectCellInfoPop

%% Estimate spot detection threshold
% define number of spots in each frame to extract score for good and bad
% spots. This is used to define a score that can be used to detect foci

win_size = 4;               % detection window size in pixels (also for bandpass filter)
min_spot_dist = 4;          % minimum distance between spots in pixels

% get cells file
[filename,pathname] = uigetfile('*.mat', 'Select','Select','MultiSelect', 'on');
if ~iscell(filename)
    name = filename;
    filename = cell(1);
    filename{1} = name;
end

figure;
good_score = [];
bad_score = [];
for hh = 1:length(filename)
    load([pathname filename{hh}])
    for ii = 1:length(cells) % cells
        for jj = 1:length(cells(ii).area) % time points
            if length(cells(ii).area) == 1
                image = cells(ii).fluo1;
            else
                image = cells(ii).fluo1{jj};
            end
            if ~isempty(image)
                % detect spots by filtering, finding local maxima and
                % fitting gaussian
                [locX, locY, fitSigma, bgsubtrInt, fitScore, score] = detectSpots(...
                    image, win_size, cells(ii).mask);

                % remove points too close to each other, remove dimmer spot
                D = squareform(pdist([locX locY]));
                % bg subtracted brightnesses
                [~,inds] = sort(bgsubtrInt,'descend');

                % loop over spots according to brighness if more than 1
                if ~isempty(D) && length(locX) > 1
                    for kk = 1:length(bgsubtrInt)
                        % other spots within threshold
                        close_inds = find(D(:,inds(kk)) < min_spot_dist & D(:,inds(kk)) ~= 0);
                        % change their values to NaN
                        if ~isempty(close_inds)
                            % change ind to NaN in different vectors
                            D(:,close_inds) = NaN;
                            D(close_inds,:) = NaN;
                            locX(close_inds) = NaN; 
                            locY(close_inds) = NaN;
                            fitSigma(close_inds) = NaN; 
                            bgsubtrInt(close_inds) = NaN; 
                            fitScore(close_inds) = NaN;
                            score(close_inds) = NaN;
                        end
                    end
                end
                
                I = image;
                size_I = size(I);
                I_crop = imcrop(I ,[win_size+1 win_size+1 size_I(2)-2.*win_size+2 size_I(1)-2.*win_size+2]);
                imagesc(I_crop)
                hold on;
                plot(locX-win_size,locY-win_size,'or')
                title([num2str(ii) ' - ' num2str(jj)])
                hold off;
                str = input('How many spots? ','s');  
                str = str2num(str);
                scores = sort(score(~isnan(score)),'descend');
                if isempty(str)
                elseif str == 0
                    bad_score = [bad_score; scores];
                elseif str < length(scores)
                    good_score = [good_score; scores(1:str)];
                    bad_score = [bad_score; scores((str+1):end)];
                else
                    good_score = [good_score; scores];
                end
            end
        end
    end
end

figure; hold on
[N, X] = hist(good_score(good_score<300),160);
N = N./sum(N);
bar(X,N)
[N, X] = hist(bad_score(bad_score>0),300);
N = N./sum(N);
bar(X,N)
xlabel('Score')

%% Spot detection by fitting a 2D gaussian
% use previous section to define the threshold

score_thres = 20;               % score treshold for a spot (Intensity * fitScore / fitSigma))
win_size = 4;                   % window size for spot region (also for bandpass filter)
mask_dilation = 2;              % amount mask dilated (in pixels) to see if spot is inside mask area
min_spot_dist = 4;              % minimum distance (in pixels) between spots to counted as separate spots

% get cells file
[filename,pathname] = uigetfile('*.mat', 'Select','Select','MultiSelect', 'on');
if ~iscell(filename)
    name = filename;
    filename = cell(1);
    filename{1} = name;
end

% figure;
for hh = 1:length(filename)
    load([pathname filename{hh}])

    for ii = 1:length(cells) % cells
        cells(ii).no_spots = NaN(1,length(cells(ii).area));
        cells(ii).spot_loc_X = [];
        cells(ii).spot_loc_Y = [];
        for jj = 1:length(cells(ii).area) % time points
            if length(cells(ii).area) == 1
                image = cells(ii).fluo1;
            else
                image = cells(ii).fluo1{jj};
            end
            if ~isempty(image)
                % detect spots by filtering, finding local maxima and
                % fitting gaussian
                [locX, locY, fitSigma, Ints, fitScore, score] = detectSpots(...
                    image, win_size, cells(ii).mask);
                
                % remove low quality spots by score
                mask_score = score > score_thres;
                locX = locX(mask_score); 
                locY = locY(mask_score);
                fitSigma = fitSigma(mask_score); 
                Ints = Ints(mask_score); 
                fitScore = fitScore(mask_score);
                score = score(mask_score);
                
                if ~isempty(locX)
                    % remove points too close to each other, remove dimmer spot
                    D = squareform(pdist([locX locY]));
                    % bg subtracted brightnesses
                    [~,inds] = sort(Ints,'descend');

                    % loop over spots according to brighness if more than 1
                    if ~isempty(D) && length(locX) > 1
                        for kk = 1:length(Ints)
                            % other spots within threshold
                            close_inds = find(D(:,inds(kk)) < min_spot_dist & D(:,inds(kk)) ~= 0);
                            % change their values to NaN
                            if ~isempty(close_inds)
                                % change ind to NaN in different vectors
                                D(:,close_inds) = NaN;
                                D(close_inds,:) = NaN;
                                locX(close_inds) = NaN; 
                                locY(close_inds) = NaN;
                                fitSigma(close_inds) = NaN; 
                                Ints(close_inds) = NaN; 
                                fitScore(close_inds) = NaN;
                                score(close_inds) = NaN;
                            end
                        end
                    end

                    % dilate masks and get boundary
                    if length(cells(ii).area) == 1
                        bw_mask = imdilate(cells(ii).mask, strel('disk', mask_dilation));
                    else
                        bw_mask = imdilate(cells(ii).mask{jj}, strel('disk', mask_dilation));
                    end
                    boundary = bwboundaries(bw_mask,8,'noholes');

                    % see if inside a cell mask
                    in = inpolygon(locX,locY,boundary{1}(:,2),boundary{1}(:,1));

                    % number of spots
                    cells(ii).no_spots(jj) = sum(in);
                    cells(ii).spot_loc_X = locX(in);
                    cells(ii).spot_loc_Y = locY(in);
                else
                    cells(ii).no_spots(jj) = 0;
                end
            end
        end
    end 

    % save with same filename
    save([pathname filename{hh}],'cells','-v7.3')
end

%% Calculate number of separate chromosome areas
% Estimates this based on HU labelled area using Otsu's thresholding.
% Limits of nucleoid area are extracted from separate exp.

min_nuc_area = 50;      % minimum area of nucleoid in pixels to be considered
max_nuc_area = 800;     % maximum area of a single nucleoid

% get cells file
[filename,pathname] = uigetfile('*.mat', 'Select','Select','MultiSelect', 'on');
if ~iscell(filename)
    name = filename;
    filename = cell(1);
    filename{1} = name;
end

for hh = 1:length(filename)   
    load([pathname filename{hh}])
	for ii = 1:length(cells) 
        % segment nucleoid using otsu
        cell2thresh = immultiply(cells(ii).fluo1,cells(ii).mask);
        thresh = multithresh(cell2thresh,2);
        mask = cell2thresh > max(thresh);
        % remove too small areas
        mask2 = bwareaopen(mask, min_nuc_area);
        % unique labels
        L = bwlabel(mask2);

        % save info to cells
        if sum(L(:)) > max_nuc_area
            cells(ii).no_chrom = NaN;
            cells(ii).nuc_mask = L;
        else
            cells(ii).no_chrom = length(unique(L))-1;
            cells(ii).nuc_mask = L;
        end
	end
    % save with same filename
    save([pathname filename{hh}],'cells','-v7.3')
end


%% Plot fluorescence concentrations

fluo1_conc = [];
fluo2_conc = [];
area = [];

% get cells file
[filename,pathname] = uigetfile('*.mat', 'Select','Select','MultiSelect', 'on');
if ~iscell(filename)
    name = filename;
    filename = cell(1);
    filename{1} = name;
end

for hh = 1:length(filename)   
    load([pathname filename{hh}])
    for ii = 1:length(cells)
        if cells(ii).alone == 1 && cells(ii).no_chrom == 1
            area = [area cells(ii).area];
            if isfield(cells,'fluo1')
                fl1 = cells(ii).fluo1 - cells(ii).bg_fl1;
                tot_fluo1 = sum(fl1(cells(ii).mask));
                fluo1_conc = [fluo1_conc tot_fluo1./cells(ii).area];
            end
            if isfield(cells,'fluo2')
                fl2 = cells(ii).fluo2 - cells(ii).bg_fl2;
                tot_fluo2 = sum(fl2(cells(ii).mask));
                fluo2_conc = [fluo2_conc tot_fluo2./cells(ii).area];
            end
        end
    end
end

figure;plot(area,fluo1_conc,'.r')
xlabel('Cell area (\mum^2)')
ylabel('Concentration (AU)')
figure;plot(area,fluo2_conc,'.b')
xlabel('Cell area (\mum^2)')
ylabel('Concentration (AU)')

xval = area;

% choose fluo1 or fluo2 as y-axis
yval = fluo1_conc;
% yval = fluo2_conc;

% add binned mean
interval = 0.5;         % interval of the mean
range = 2:interval:10;  % range of the mean

means = zeros(1,length(range)-1);
sems = zeros(1,length(range)-1);
for ii = 2:length(range)
    mask = xval >= range(ii-1) & xval < range(ii);
    means(ii-1) = mean(yval(mask));
    sems(ii-1) = std(yval(mask))./sqrt(length(yval(mask)));
end
hold on;errorbar(range(1:(end-1))+interval/2,means,sems,'k','LineWidth',1)

hold on;plot(mean(area),mean(fluo1_conc),'ok')
hold on;errorbar(mean(area),mean(fluo1_conc),std(fluo1_conc)./sqrt(length(fluo1_conc)))

% plot KS density
figure; hold on;
density = ksdensity([area',fluo1_conc'], [area',fluo1_conc']);
scatter(area', fluo1_conc', 20, density,'filled');
xlabel('Cell area (\mum^2)')
ylabel('Total HupA-mCherry intensity (AU)')

[f, xi, bw] = ksdensity([area',fluo1_conc'], [area',fluo1_conc']);
figure;ksdensity([area',fluo1_conc'], [area',fluo1_conc'],'Bandwidth',bw,'PlotFcn','contour');

figure; hold on;
density = ksdensity([fluo1_conc',fluo2_conc'], [fluo1_conc',fluo2_conc']);
scatter(fluo1_conc', fluo2_conc', 20, density,'filled');
xlabel('HupA-mCherry concentration (AU)')
ylabel('RpsB concentration (AU)')

figure; hold on;
density = ksdensity([fluo1_conc',fluo2_conc'], [fluo1_conc',fluo2_conc']);
scatter(fluo1_conc', fluo2_conc', 20, cell_area,'filled');
xlabel('HupA-mCherry concentration (AU)')
ylabel('RpsB concentration (AU)')

figure;hist(cell_area,100)

%% Plotting fluorescence concentrations with repeats

number_repeats = 3;         % number of independent experiments

limits = [2.2 10];          % plotting limits
interval = 0.5;             % interval for binning data
range = (limits(1)-interval/2):interval:(limits(2)+interval*1.5);    % range for binning data
area_bins = range(1:(end-1))+interval/2;

means_all = NaN(number_repeats, length(range)-1);
fluo1_conc_all = [];
fluo2_conc_all = [];
area_all = [];
for ss = 1:number_repeats
    fluo1_conc = [];
    fluo2_conc = [];
    area = [];

    % get cells file
    [filename,pathname] = uigetfile('*.mat', 'Select','Select','MultiSelect', 'on');
    if ~iscell(filename)
        name = filename;
        filename = cell(1);
        filename{1} = name;
    end
    for hh = 1:length(filename)   
        load([pathname filename{hh}])
        for ii = 1:length(cells)
            if cells(ii).alone == 1 && cells(ii).no_chrom == 1
                area = [area cells(ii).area];
                fl1 = cells(ii).fluo1 - cells(ii).bg_fl1;
                tot_fluo1 = sum(fl1(cells(ii).mask));
                fluo1_conc = [fluo1_conc tot_fluo1./cells(ii).area];
                if isfield(cells,'fluo2')
                    fl2 = cells(ii).fluo2 - cells(ii).bg_fl2;
                    tot_fluo2 = sum(fl2(cells(ii).mask));
                    fluo2_conc = [fluo2_conc tot_fluo2./cells(ii).area];
                end
            end
        end
    end
    % choose which fluorescence channel
    xval = area;
    yval = fluo1_conc;
%     yval = fluo2_conc;

    % add moving mean
    means = zeros(1,length(range)-1);
    sems = zeros(1,length(range)-1);
    for ii = 2:length(range)
        mask = xval >= range(ii-1) & xval < range(ii);
        means(ii-1) = mean(yval(mask));
        sems(ii-1) = std(yval(mask))./sqrt(length(yval(mask)));
    end
    
    means_all(ss,:) = means;
    fluo1_conc_all = [fluo1_conc_all fluo1_conc];
    fluo2_conc_all = [fluo2_conc_all fluo2_conc];
    area_all = [area_all area];
    disp(['number of cells ' num2str(length(yval))])
end

figure;plot(area_all,fluo1_conc_all,'.b')
xlabel('Cell area (\mum^2)')
ylabel('Concentration (AU)')
figure;plot(area_all,fluo2_conc_all,'.b')
xlabel('Cell area (\mum^2)')
ylabel('Concentration (AU)')

means = mean(means_all,'omitnan');
stds = std(means_all,'omitnan');
sems = std(means_all,'omitnan')./sqrt(number_repeats);

hold on;plot(area_bins,means,'k','LineWidth',1)
hold on;plot(area_bins,means-stds,'k','LineWidth',0.5)
hold on;plot(area_bins,means+stds,'k','LineWidth',0.5)

%check quantiles
limits = quantile(area_all,[0.025 0.975]);
    
%% Plotting total intensity and concentration vs cell area

% join different XY cell together
[filename,pathname] = uigetfile('*.mat', 'Select','MultiSelect', 'on');

cell_areas = [];
sum_ints = [];
n = 1;

if ~iscell(filename)
    filename2 = filename;
    filename = cell(1,1);
    filename{1} = filename2;
end

for kk = 1:length(filename)
    load([pathname filename{kk}]);
    cell_areas = [cell_areas [cells.area]];
    sum_ints = [sum_ints [cells.fluo1_int]];
end
conc_ints = sum_ints./cell_areas;

figure; hold on;
density = ksdensity([cell_areas',sum_ints'], [cell_areas',sum_ints']);
scatter(cell_areas', sum_ints', 20, density,'filled');
xlabel('Cell area (\mum^2)')
ylabel('Total HupA-mCherry intensity (AU)')

figure; hold on;
density = ksdensity([cell_areas',conc_ints'], [cell_areas',conc_ints']);
scatter(cell_areas', conc_ints', 20, density,'filled');
xlabel('Cell area (\mum^2)')
ylabel('HupA-mCherry concentration (AU)')

% fit int to a line
p = polyfit(cell_areas,sum_ints,3);
x1 = linspace(min(cell_areas),max(cell_areas));
y1 = polyval(p,x1);
hold on
plot(x1,y1)

% fit concentration to a curve
p = polyfit(cell_areas(cell_areas > 10),conc_ints(cell_areas > 10),1);
x1 = linspace(min(cell_areas),30);
y1 = polyval(p,x1);
hold on
plot(x1,y1)
