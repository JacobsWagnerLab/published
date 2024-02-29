%% HU nucleoid detection accuracy compared to ParB spot detection
% assumes fluo1 is ParB spots and fluo2 is HU fluorescence

% collect info
collectCellInfoPop

%% Spot detection with ParB
score_thres = 50; %20              % score treshold for a spot (Intensity * fitScore / fitSigma))
win_size = 4;                   % window size for spot region (also for bandpass filter)
mask_dilation = 1;              % amount mask dilated (in pixels) to see if spot is inside mask area
min_spot_dist = 5;              % minimum distance (in pixels) between spots to counted as separate spots

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
        cells(ii).spots = struct;
        for jj = 1:length(cells(ii).area) % time points
            cells(ii).spots.locX = [];
            cells(ii).spots.locY = [];
            cells(ii).spots.fitSigma = [];
            cells(ii).spots.Int = [];
            cells(ii).spots.fitScore = [];
            cells(ii).spots.score = [];
            cells(ii).no_spots = [];
            
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
                    cells(ii).spots.locX = locX(in);
                    cells(ii).spots.locY = locY(in);
                    cells(ii).spots.fitSigma = fitSigma(in);
                    cells(ii).spots.Int = Ints(in);
                    cells(ii).spots.fitScore = fitScore(in);
                    cells(ii).spots.score = score(in);
                    cells(ii).no_spots = length(locX(in));
                else
                    cells(ii).no_spots = 0;
                end
            end
        end
    end 

    % save with same filename
    save([pathname filename{hh}],'cells')
end

%% Detect number of nucleoids with HU

min_cell_area = 7;      % minimum area of a cell
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
        cell2thresh = immultiply(cells(ii).fluo2,cells(ii).mask);
        thresh = multithresh(cell2thresh,2);
        mask = cell2thresh > max(thresh);
        mask2 = bwareaopen(mask, min_nuc_area);
        L = bwlabel(mask2);
        
        % save info to cells
        if sum(L(:)) > max_nuc_area || cells(ii).area < min_cell_area
            cells(ii).no_nucleoid = NaN;
            cells(ii).nucleoid_mask = L;
        else
            cells(ii).no_nucleoid = length(unique(L))-1;
            cells(ii).nucleoid_mask = L;
        end
    end
    % save with same filename
    save([pathname filename{hh}],'cells')
end

%% Check accuracy

ParB_no = [];
for ii = 1:length(cells)
    if cells(ii).no_nucleoid == 1
        if ~isnan(cells(ii).no_spots) && cells(ii).no_spots ~= 0
            ParB_no = [ParB_no cells(ii).no_spots];
        end
    end
end
sum(ParB_no == 1)./length(ParB_no)

% Correct wrong ones
figure; 
for ii = 1:length(cells)
    if (cells(ii).no_nucleoid == 1 && cells(ii).no_spots ~= 1)
        tiledlayout(1,3)
        nexttile
        imagesc(cells(ii).fluo1) 
        title(num2str(cells(ii).spots.score))
        nexttile
        imagesc(cells(ii).fluo2)
        title(num2str(ii))
        nexttile
        imagesc(cells(ii).nucleoid_mask)
        str = input('Number of ParB spots: ','s');  
        str = str2num(str);
        if ~isempty(str)
            cells(ii).no_spots = str; 
        end
    end
end

%% Single chromosome area
areas = [];
for ii = 1:length(cells)
    if cells(ii).no_spots == 2
        areas = [areas sum(cells(ii).nucleoid_mask(:))];
    end
end

figure; hist(areas,100)
