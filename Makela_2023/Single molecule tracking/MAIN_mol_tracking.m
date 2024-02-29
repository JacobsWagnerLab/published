%% Single molecule diffusion analysis
% use the following pipeline to obtain tracks for single molecules

%%% 1. Detect cell outlines using MicrobeTracker:
microbeTracker

%%% 2. Localize single molecules in each frame:
open detect_mol

%%% 3. Track molecules inside cell regions:
open track_mol

% 7/7/22 Jarno Makela

%% Parameters for experiments
% used for calculations below
params = struct;
params.pixel = 0.16;                    % length per pixel in um
params.dT = 0.0107;                     % time interval in seconds
params.DhistMinSteps = 10;              % minimum number of steps for a track
params.truncated_tracks = 0;            % (not used) truncate tracks to the same length 
params.sigmaNoise = 0;                  % (not used) localization noise 

%% Measure D for each molecule

% get *filt_locoli.mat files
[filename,pathname] = uigetfile('*.mat', 'Select','Select','MultiSelect', 'on');
if ~iscell(filename)
    name = filename;
    filename = cell(1);
    filename{1} = name;
end

D = cell(1,10000);      % pre-allocate variable
no_cells = 0;
ind = 1;
for hh = 1:length(filename)   
    load([pathname filename{hh}])
    no_cells = no_cells + size(movie_data.cellROI_data,1);

    for ii = 1:size(movie_data.cellROI_data,1)
        tracks = movie_data.cellROI_data(ii,1).tracks;
        if isempty(tracks) == 0
            % use histDv2 function to calculate D values
            [D_single, long_enough_tracks] = histDv2(tracks,params);
        end
        % save data into cell
        D{ind} = D_single;
        ind = ind + 1;
    end
    disp('.')
end
% convert cell into a vector
D = cat(1, D{:});

disp(['Number of cells = ' num2str(no_cells)])
disp(['Mean D* value = ' num2str(mean(D))])
disp(['Number of molecules = ' num2str(length(D))])

% plot distribution
no_bins = 60;
[~,edges] = histcounts(log10(D),no_bins);
[N, X] = hist(D,10.^edges);
N = N./((sum(N)*(log10(X(2))-log10(X(1))))); 

figure;
bar(log10(X),N,'BarWidth',1,'FaceColor',[.9 .9 .9]);
xlabel('log_{10}(D*)(\mum^2/s)')
ylabel('Probability density')
axis([-2.5 1 0 Inf])
hold on;
plot(log10(X),N,'LineWidth', 1)

% ksdensity of the data
[f,xi] = ksdensity(log10(D)); 
figure;plot(xi,f,'k','linewidth', 2);
axis([-2.5 1 0 Inf])

%% EM fit of gaussian mixture model
% Expectation maximization algorithm for fitting model to data
% Note! Initiation of search is random and can give different local maxima

% run previous section before this

positions = [NaN NaN 0.4]; 	% number of fitted species; numbers indicate fixed locations of peaks, NaN is used for free location
x = -2.5:0.01:1;            % range of values plotted for GMM fit
low_perc = 0.005;           % percentage of smallest values ignored in fitting

% convert D values into logarithmic space for fitting
values = log10(D');         

% remove smallest X% for better fitting
lower = quantile(values, low_perc);
values = values(values > lower);

% fit GMM to data
[mu_est, sigma_est, w_est, counter, difference] = gmm_fixed_param(values, length(positions), positions, 1.0e-3);

% plot estimated distributions
figure;
hold on;
p_est_all = [];
for ii = 1:length(mu_est)
    % multiply normal density dist with weight
    p_est_single = w_est(ii) * norm_density(x, mu_est(ii), sigma_est(ii));
    plot(x, p_est_single, 'k--', 'linewidth', 2);
    % keep densities
    p_est_all = [p_est_all; p_est_single];
    disp(['Mean value ' num2str(mu_est(ii)) ' - weight ' num2str(w_est(ii))])
end
% sort GMM classes from slowest diffusion to largest (class 1 bound)
[mu_est,I] = sort(mu_est);
sigma_est = sigma_est(I);
w_est = w_est(I);

plot(x, sum(p_est_all), 'k-', 'linewidth', 2);
title(['Mean value ' num2str(mu_est) ', weight ' num2str(w_est')])

% save the GMM
save('GMM_model.mat','mu_est','sigma_est','w_est')

%% Estimate bound fraction in individual cells based on GMM
% load GMM with mu_est, sigma_est, w_est (from previous section)

min_no_mol = 50;        % minimum number of molecules per cell

areas = [];             % cell areas
fracs = [];             % fraction of molecules bound in a cell
mol_no = [];            % total number of molecules
indeces = [];           % cell indeces
    
% get *filt_locoli.mat files
[filename,pathname] = uigetfile('*.mat', 'Select','Select','MultiSelect', 'on');
if ~iscell(filename)
    name = filename;
    filename = cell(1);
    filename{1} = name;
end

for hh = 1:length(filename)   
    load([pathname filename{hh}])

    % initialize new fields
    movie_data.cellROI_data(1).D = [];
    movie_data.cellROI_data(1).frac_bound = [];
    for ii = 1:size(movie_data.cellROI_data,1)
        tracks = movie_data.cellROI_data(ii,1).tracks;
        if isempty(tracks) == 0
            % use histDv2 function to calculate D values
            [D, long_enough_tracks] = histDv2(tracks,params);
        end
        movie_data.cellROI_data(ii).D = D;
        
        % sort GMM classes from smallest diffusion to largest (class 1 is bound)
        [mu_est,I] = sort(mu_est);
        sigma_est = sigma_est(I);
        w_est = w_est(I);
        
        % classify log10 D values according to GMM probability
        values = log10(D'); 
        class = zeros(length(mu_est),length(values));
        for jj = 1:length(mu_est)
            class(jj, :) = w_est(jj) * norm_density(values, mu_est(jj), sigma_est(jj));
        end

        % normalize data by columns
        class = class ./ repmat(sum(class), length(mu_est), 1);
        % largest probability defines the class of a track
        [~,mol_species] = max(class,[],1);
        
        % find out fraction in class 1 that corresponds to bound molecules
        movie_data.cellROI_data(ii).frac_bound = sum(mol_species == 1)./length(mol_species);
        
        % collect data
        if ~isempty(movie_data.cellROI_data(ii).D) && ...
                length(movie_data.cellROI_data(ii).D) >= min_no_mol
            areas = [areas movie_data.cellROI_data(ii).area.*params.pixel^2];
            fracs = [fracs movie_data.cellROI_data(ii).frac_bound];
            mol_no = [mol_no length(movie_data.cellROI_data(ii).D)];
            indeces = [indeces ii];
        end
    end
end

% plotting
figure; 
hold on;
plot(areas, fracs,'.','Color',[0.8 0.8 0.8],'MarkerSize',12)
axis([0 10 0 1])
xlabel('Cell area (\mum^2)')
ylabel('Fraction of RNAP bound')
title(['No. cells ' num2str(length(areas))])

figure;scatter(areas, fracs,[],mol_no)

hold on;errorbar(mean(areas),mean(fracs),std(fracs)./sqrt(length(fracs)))

% 95% confidence interval from bootstrapping
centers = [1.5 3];              % range for binning
interval = 0.2;                 % interval for binning data

range = (centers(1)-interval/2):interval:(centers(2)+interval/2);    % range for binning data
area_bins = range(1:(end-1))+interval/2;

% find bins
[~,~,idx] = histcounts(areas,range);
upper = NaN(1,length(1:max(idx)));
lower = NaN(1,length(1:max(idx)));
mean_CI = NaN(1,length(1:max(idx)));
for ii = 1:max(idx)
    % bootstrap 10k times
    m = bootstrp(10000,@mean,fracs(idx==ii));
    y = quantile(m,[0.025 0.975]);
    upper(ii) = y(2);
    lower(ii) = y(1);
    mean_CI(ii) = mean(m);
end

% plot 95% CI
hold on;plot(area_bins,mean_CI,'k','LineWidth',1)
hold on;plot(area_bins,upper,'k','LineWidth',0.5)
hold on;plot(area_bins,lower,'k','LineWidth',0.5)

%% Ksdensity of D for cells with different areas
% bin cells according to their cell area

bin_centers = 2:1:10;           % bin centers
bin_width = 0.5;                % width for one direction (total 2*width)

% get *filt_locoli.mat files
[filename,pathname] = uigetfile('*.mat', 'Select','Select','MultiSelect', 'on');
if ~iscell(filename)
    name = filename;
    filename = cell(1);
    filename{1} = name;
end

for hh = 1:length(filename)   
    load([pathname filename{hh}])
    % initialize D field
    movie_data.cellROI_data(1).D = [];
    movie_data.cellROI_data(1).file = [];
    for ii = 1:size(movie_data.cellROI_data,1)
        tracks = movie_data.cellROI_data(ii,1).tracks;
        D = [];
        if isempty(tracks) == 0
            % use histDv2 function to calculate D values
            [D, long_enough_tracks] = histDv2(tracks,params);
        end
        movie_data.cellROI_data(ii).D = D';
        movie_data.cellROI_data(ii).file = hh;
    end
    % collect data from multiple files
    if hh == 1
        data = movie_data.cellROI_data;
    else
        data = [data; movie_data.cellROI_data];
    end  
end

% find data indeces with correct cell areas
areas_all = [data.area].*params.pixel^2;
bin_indeces = cell(1,length(bin_centers));
logD_bins = cell(1,length(bin_centers));
ksd_f = cell(1,length(bin_centers));
ksd_xi = cell(1,length(bin_centers));
for ii = 1:length(bin_centers)
    bin_indeces{ii} = find((areas_all >= bin_centers(ii) - bin_width) & (areas_all <= bin_centers(ii) + bin_width)); 
    values = log10([data(bin_indeces{ii}).D]);         % convert D values into logarithmic space
    logD_bins{ii} = values;
    [ksd_f{ii},ksd_xi{ii}] = ksdensity(values);
end

% plot ksdensity
figure; hold on;
cm = colormap(parula(length(ksd_xi)));  % defines colormap for X number of elements
for ii = 1:length(ksd_f)
    plot(ksd_xi{ii},ksd_f{ii},'Color',cm(ii,:), 'LineWidth',2);
end
axis([-2.5 1 0 Inf])

%% EM fit of GMM model based on cell sizes binning
% Only include cells in each bin for EM fitting
% Bin centers and width can be independently tuned

bin_centers = 2.5:1:10;         % bin centers
bin_width = 0.5;                % width of bin from center to both directions
positions = [NaN NaN 0.42];     % fixed locations of peaks, NaN is not defined
x = -2.5:0.01:1;                % range of values plotted for mixture model
low_perc = 0.005;               % percentage of lower part ignored in fitting

% get *filt_locoli.mat files
[filename,pathname] = uigetfile('*.mat', 'Select','Select','MultiSelect', 'on');
if ~iscell(filename)
    name = filename;
    filename = cell(1);
    filename{1} = name;
end

% collect all data
for hh = 1:length(filename)   
    load([pathname filename{hh}])
    % initialize D field
    movie_data.cellROI_data(1).D = [];
    for ii = 1:size(movie_data.cellROI_data,1)
        tracks = movie_data.cellROI_data(ii,1).tracks;
        D = [];
        if isempty(tracks) == 0
            % use histDv2 function to calculate D values
            [D, long_enough_tracks] = histDv2(tracks,params);
        end
        movie_data.cellROI_data(ii).D = D';
    end
    % collect data from multiple files
    if hh == 1
        data = movie_data.cellROI_data;
    else
        data = [data; movie_data.cellROI_data];
    end  
end
% cell areas
areas = [data.area].*params.pixel^2;

% find data indeces with correct cell areas and estimate GMM for all bins
bin_indeces = cell(1,length(bin_centers));
mu_est = cell(1,length(bin_centers));
sigma_est = cell(1,length(bin_centers));
w_est = cell(1,length(bin_centers));
logD_bins = cell(1,length(bin_centers));
for ii = 1:length(bin_centers)
    bin_indeces{ii} = find((areas >= bin_centers(ii) - bin_width) & (areas <= bin_centers(ii) + bin_width)); 
    values = log10([data(bin_indeces{ii}).D]);         % convert D values into logarithmic space for fitting
    logD_bins{ii} = values;
    
    if ~isempty(values)
        % remove smallest X% for better fitting
        lower = quantile(values, low_perc);
        values = values(values > lower);

        % fit GMM
        [mu_est{ii}, sigma_est{ii}, w_est{ii}, counter, difference] = ...
            gmm_fixed_param(values, length(positions), positions, 1.0e-3);

        % sort GMM classes from smallest diffusion to largest (class 1 bound)
        [mu_est{ii},I] = sort(mu_est{ii});
        sigma_est{ii} = sigma_est{ii}(I);
        w_est{ii} = w_est{ii}(I);
    end
end

% plot D* and GMM of a specific bin
bin_no = 6;

% plot histogram
figure; hold on;
[~,edges] = histcounts(logD_bins{bin_no},60);
[N, X] = hist(10.^logD_bins{bin_no},10.^edges);
N = N./((sum(N)*(log10(X(2))-log10(X(1))))); 
bar(log10(X),N,'BarWidth',1,'FaceColor',[.9 .9 .9]);
xlabel('log_{10}(D*)(\mum^2/s)')
ylabel('Probability density')

% plot GMM model
p_est_all = [];
for ii = 1:length(mu_est{bin_no})
    p_est_single = w_est{bin_no}(ii) * norm_density(x, mu_est{bin_no}(ii), sigma_est{bin_no}(ii));
    plot(x, p_est_single, 'k--', 'linewidth', 2);
    p_est_all = [p_est_all; p_est_single];
end
plot(x, sum(p_est_all), 'k-', 'linewidth', 2);
title(['Bin no. ' num2str(bin_no) ', mean value ' num2str(mu_est{bin_no}(1)) ', weight ' ...
    num2str(w_est{bin_no}(1)) ' n=' num2str(length(logD_bins{bin_no}))])

% save GMM for the future use
save('GMM_model.mat','mu_est', 'sigma_est', 'w_est', 'bin_centers', 'bin_width', 'low_perc')

%% Estimate bound fraction in individual cells based on GMM with size bins
% load GMM model with mu_est, sigma_est and w_est

min_no_mol = 50;            % minimum number of tracks per cell

% params for binning mean
borders = [2 10];           % bin locations
interval = 0.5;             % interval for binning data
range = (borders(1)-interval/2):interval:(borders(end)+interval/2);    % range for binning data
area_bins = range(1:(end-1))+interval/2;

% get *filt_locoli.mat files
[filename,pathname] = uigetfile('*.mat', 'Select','Select','MultiSelect', 'on');
if ~iscell(filename)
    name = filename;
    filename = cell(1);
    filename{1} = name;
end

for hh = 1:length(filename)   
    load([pathname filename{hh}])
    % initialize D field
    movie_data.cellROI_data(1).D = [];
    movie_data.cellROI_data(1).file = [];
    for ii = 1:size(movie_data.cellROI_data,1)
        tracks = movie_data.cellROI_data(ii,1).tracks;
        D = [];
        if isempty(tracks) == 0
            % use histDv2 function to calculate D values
            [D, long_enough_tracks] = histDv2(tracks,params);
        end
        movie_data.cellROI_data(ii).D = D';
        movie_data.cellROI_data(ii).file = hh;
    end
    % collect data from multiple files
    if hh == 1
        data = movie_data.cellROI_data;
    else
        data = [data; movie_data.cellROI_data];
    end  
end

% classify individual tracks according to cell size specific GMM
areas = [];
fracs = [];
indeces = [];
mol_no = [];
data(1).frac_bound = [];
for ii = 1:length(data)
    % cells that have enough tracks
    if ~isempty(data(ii).D) && length(data(ii).D) >= min_no_mol  
        % check if GMM is defined for each bin separately
        if ~iscell(mu_est)
            % classify log10 D values according to GMM probability
            values = log10(data(ii).D'); 
            class = zeros(length(mu_est),length(values));
            for jj = 1:length(mu_est)
                class(jj,:) = w_est(jj) * norm_density(values, mu_est(jj), sigma_est(jj));
            end

            % normalize data by columns
            class = class ./ repmat(sum(class), length(mu_est), 1);
            % largest probability defines the class
            [~,mol_species] = max(class,[],1);
        else
            % find closest bin center
            [~, ind] = min(abs(bin_centers-(data(ii).area.*params.pixel^2)));
            
            % classify log10 D values according to GMM probability
            values = log10(data(ii).D'); 
            class = zeros(length(mu_est{ind}),length(values));
            for jj = 1:length(mu_est{ind})
                class(jj,:) = w_est{ind}(jj) * norm_density(values, mu_est{ind}(jj), sigma_est{ind}(jj));
            end

            % normalize data by columns
            class = class ./ repmat(sum(class), length(mu_est{ind}), 1);
            % largest probability defines the class
            [~,mol_species] = max(class,[],1);
        end

        % find out fraction in class 1 that corresponds to bound molecules
        data(ii).frac_bound = sum(mol_species == 1)./length(mol_species);
        mol_no = [mol_no length(data(ii).D)];
        areas = [areas data(ii).area.*params.pixel^2];
        fracs = [fracs data(ii).frac_bound];
        indeces = [indeces ii];
    end 
end

% plotting
figure; 
plot(areas, fracs,'.','Color',[0.8 0.8 0.8],'MarkerSize',12)
xlabel('Cell area (\mum^2)')
ylabel('Fraction of RNAP bound')
axis([0 10 0 1])
title(['No. cells ' num2str(length(areas))])

figure; 
scatter(areas, fracs, [], mol_no)
xlabel('Cell area (\mum^2)')
ylabel('Fraction of RNAP bound')

% bin step_size data and calculate mean, std, sem for fraction plot
[~,~,idx] = histcounts(areas,range);
means = accumarray(idx(idx > 0)',fracs(idx > 0),[],@mean);
stds = accumarray(idx(idx > 0)',fracs(idx > 0),[],@std);
sems = accumarray(idx(idx > 0)',fracs(idx > 0),[],@(x) std(x)./sqrt(length(x)));

% plot SEM
hold on;plot(area_bins,means,'k','LineWidth',1)
hold on;plot(area_bins,means-sems,'k','LineWidth',0.5)
hold on;plot(area_bins,means+sems,'k','LineWidth',0.5)

% 95% confidence interval from bootstrap
[~,~,idx] = histcounts(areas,range);
upper = NaN(1,length(1:max(idx)));
lower = NaN(1,length(1:max(idx)));
mean_CI = NaN(1,length(1:max(idx)));
for ii = 1:max(idx)
    m = bootstrp(10000,@mean,fracs(idx==ii));
    y = quantile(m,[0.025 0.975]);
    upper(ii) = y(2);
    lower(ii) = y(1);
    mean_CI(ii) = mean(m);
end

% plot 95% CI
hold on;plot(area_bins,mean_CI,'k','LineWidth',1)
hold on;plot(area_bins,upper,'k','LineWidth',0.5)
hold on;plot(area_bins,lower,'k','LineWidth',0.5)

%% Estimate the number of RNAP bound to chromosome
% run the previous section before

% RNAP intensity data
fluo1_int = [];
fluo2_int = [];
fl_area = [];

% get files for total fluorescence data
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
            fl_area = [fl_area cells(ii).area];
            fl1 = cells(ii).fluo1 - cells(ii).bg_fl1;
            fluo1_int = [fluo1_int sum(fl1(cells(ii).mask))];
            if isfield(cells,'fluo2')
                fl2 = cells(ii).fluo2 - cells(ii).bg_fl2;
                fluo2_int = [fluo2_int sum(fl2(cells(ii).mask))];
            end
        end
    end
end
% in case fluo2 does not exist
if isempty(fluo2_int)
    fluo2_int = fluo1_int;
end

figure;plot(fl_area,fluo2_int,'.b')
xlabel('Cell area (\mum^2)')
ylabel('Total intensity (AU)')

% bin the data and estimate 95% confidence interval for mean
% 95% confidence interval from bootstrap
[~,~,idx] = histcounts(fl_area,range);
fl_upper = NaN(1,length(1:max(idx)));
fl_lower = NaN(1,length(1:max(idx)));
fl_mean_CI = NaN(1,length(1:max(idx)));
for ii = 1:max(idx)
    % bootstrap 10k times
    m = bootstrp(10000,@mean,fluo2_int(idx==ii));
    y = quantile(m,[0.025 0.975]);
    fl_upper(ii) = y(2);
    fl_lower(ii) = y(1);
    fl_mean_CI(ii) = mean(m);
end

hold on;plot(area_bins,fl_mean_CI,'r','LineWidth',1)
hold on;plot(area_bins,fl_upper,'r','LineWidth',0.5)
hold on;plot(area_bins,fl_lower,'r','LineWidth',0.5)

% estimate bound RNAP quantity from total int and bound fraction
RNAP_bound = fl_mean_CI.*mean_CI;

% percentage of 95% CI of the mean
err_frac = (upper-lower)./mean_CI;
err_fl = (fl_upper-fl_lower)./fl_mean_CI;
% combined error estimate converted to actual values
err = RNAP_bound.*sqrt(err_frac.^2+err_fl.^2);

figure;plot(area_bins,RNAP_bound,'r','LineWidth',1)
hold on;plot(area_bins,RNAP_bound+(err/2),'r','LineWidth',0.5)
hold on;plot(area_bins,RNAP_bound-(err/2),'r','LineWidth',0.5)
