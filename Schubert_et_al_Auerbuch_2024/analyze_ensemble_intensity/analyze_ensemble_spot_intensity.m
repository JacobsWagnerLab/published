%% Analysis of Yersinia psueotuberculosis plasmid localization data. 
% Made by Joshua W. McCausland in Christine Jacobs-Wagner's lab, August
% 2024.
% This script will detail how to pool all of Jessica's imaging data,
% analyze the ensemble intensity, and quantify the number of plasmid foci
% in cells. The script is split into several sections that have separate
% functions. 
%   - Section 1. Prepare the main dataframe. I run this first. If running
%   from Github, then you can skip this section. Future sections load the
%   dataframe that is generated here. 
%
%   - Section 2. Generate Otsu thresholds of each condition. This will
%   calculate a single Otsu threshold of dpcnB data per condition 
%   ('induced' vs. 'uninduced') and apply that to all data within said
%   condition. 
%
%   - Section 3. Comparative histograms of wt vs. pcnB using the histograms
%   of Section 2 as a base. This is for the paper. 
%
%   - Section 4. Take the Otsu thresholds of the main dataframe generated
%   in Section 1 and add a new field: whether or not any particular cell is
%   below the Otsu threshold or not. Update and save the main dataframe.
%
%   - Section 5. Plot a histogram for each replicate of each condition.
%   This is for visual inspection to see how each replicate varies, and
%   where the calculated Otsu threshold falls.
%
%   - Section 6. Calculate the fraction of cell that are GFP-negative. This
%   uses the fraction of cells below the threshold for each replicate and
%   plots them as a bar chat, each dot representing a replicate.
%
%   - Section 7. The fraction of negative cells that are GFP negative for
%   each collective condition. This does not consider a per-replicate case,
%   rather the ensemble results based on the otsu results of Sections 2/3. 
%
%   - Section 8. Binomial distribution plots of the number of spots per
%   condition. Cells with 0 spots are excluded for this section.
%
%   - Section 9. Binomial distributions on a per-replicate basis just for
%   visual inspection. cells with 0 spots are included for this section.
%
%   - Section 10. Correlation of the fraction of cells below the otsu
%   threshold vs. the number of cells with 0 spots. 
%
%   - Section 11. Spot intensity histograms for cells whose fluorescence
%   intensity fell below the Otsu threshold.
%
%   - Section 12. Cell length violin plots comparing below/above the otsu
%   threshold for each condition.
%
%   - Section 13. Cell width violin plots comparing below/above the otsu
%   threshold for each condition.
%
%%%%%%%%%%%% REQUIRED SUPPLEMENTAL SCRIPTS FROM MATHWORKS CENTRAL %%%%%%%
% I use some supplemental scripts from Mathworks to make this run for
% dataframe construction and plotting. I reccoment you download and put
% these in  your path for this script to work best.
%   - Violin Plot version 1.7.0.0 by Holger Hoffmann. Used for sections 12
%   and 13.
%       https://www.mathworks.com/matlabcentral/fileexchange/45134-violin-plot
%
%   - Shaded Error Bar Version 1.65.0.4 by Rob Campbell. Used for section
%   3.
%       https://www.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar
%
%   - Padcat Version 1.4.0.0 by Jos. This is used for concatenating uneven
%   vectors into a matrix, useful for sections 3, 12, 13.
%       https://www.mathworks.com/matlabcentral/fileexchange/22909-padcat
%% Section 1. Prepare the main dataframe.
% Jessica has cells segmented by Oufti. Here, we will take the Oufti
% segmentations to 1) count the number of plasmid foci within each Oufti
% mesh and 2) quantify the corresponding fluorescence of these cells. Spots
% were identified via ThunderSTORM, so we can verify ThunderSTORM spot
% detection accuracy by checking detected fluorescence vs. number of
% observed spots in addition to our quantification. 
clear; clc; close all;

% Making polygons from Oufti meshes generates a benign warning message per
% cell that doesn't mean anything for me. That slows down the code as it 
% iterates through tens of thousands of cells. I turn off warning messages 
% for this chunk and turn them back on at the end. 
warning('off','all');

px_size = 0.065; %microns/px
cell_length_minimum = 1; %minimum cell length for inclusion in the dataset.

% Jessica has organized all her files into one directory for me to easily
% access. 
experiment_directory = '/Volumes/Data_02/Jessica Zhang/forJosh_ypstb';

% Identify directories for each condition ("wt_inducing",
% "wt_noninducting", "pcnb_inducing", "pcnb_noninducing")
condition_directories = dir(experiment_directory);
condition_directories = condition_directories(~contains({condition_directories(:).name},'.'));

% This counter is for iteratively appending to the intensity dataframe.
counter = 1;
for condition_index = 1:length(condition_directories)
    current_condition = condition_directories(condition_index).name;
        
        % pull up the current condition directory. 
        current_condition_directory = condition_directories(condition_index);
        % identify the oufti mesh files.
        yps_file_list = dir([current_condition_directory.folder '/' current_condition_directory.name '/*.mat']);
        % find the directories with ThunderSTORM outputs.
        ts_directory_list = dir([current_condition_directory.folder '/' current_condition_directory.name '/*ts_results*']);
        for yps_file_list_ind = 1:length(yps_file_list)
    
            % Pull up the current ThunderSTORM outputs folder.
            current_thunderSTORM_directory = ts_directory_list(yps_file_list_ind);
            % Identify all the ThunderSTORM files (one file per region)
            ts_files = dir([current_thunderSTORM_directory.folder '/' current_thunderSTORM_directory.name '/*.csv']);
            % identify the current Oufti mesh file.
            current_yps_file = yps_file_list(yps_file_list_ind);

            %grab the labeling condition from the Oufti filename.
            name_components = split(current_yps_file.name,'_');
            replicate = name_components{1};
            labeling_condition = [name_components{2} '_' name_components{3}];
            
            %Iterate through every cell and append it to a master
            %dataframe.
            oufti_file_curr = load([current_yps_file.folder '/' current_yps_file.name]);
            for idxa = 1:length(oufti_file_curr.cellList.meshData)
                region = oufti_file_curr.cellList.meshData{idxa};
                % find the ThunderSTORM file corresponding to the current
                % region.
                ts_file_curr = readtable([ts_files(idxa).folder '/' ts_files(idxa).name]);
                for idxb = 1:length(region)
                    cell_curr = region{idxb};
                    if length(cell_curr.signal1)*px_size > cell_length_minimum
                        % Calculate cell width by measuring the euclidean
                        % distance along every  point of the cell length. Then
                        % calculate the mean of the top 1/3 of these measurements.
                        width = sort(sqrt((cell_curr.mesh(:,1)-cell_curr.mesh(:,3)).^2+(cell_curr.mesh(:,2)-cell_curr.mesh(:,4)).^2),'descend');
                        width = mean(width(1:floor(length(width)/3)))*px_size;
    
                        % Construct a polygon of the current Oufti mesh
                        pgon = polyshape([cell_curr.mesh(:,1);cell_curr.mesh(:,3)]*px_size,[cell_curr.mesh(:,2);cell_curr.mesh(:,4)]*px_size);
                        
                        % Count the number of plasmid foci in the polygon using
                        % ThunderSTORM localization coordinates.
                        points_in_cell = inpolygon(ts_file_curr.x_nm_/1000,ts_file_curr.y_nm_/1000,pgon.Vertices(:,1),pgon.Vertices(:,2));
                        
                        % Add all current labels and measurements of current
                        % cell to the dataframe.
                        yps_data(counter).filename = current_yps_file.name(1:end-4);
                        yps_data(counter).region = idxa;
                        yps_data(counter).cell = idxb;
                        yps_data(counter).condition = current_condition;
                        yps_data(counter).replicate = str2num(replicate(end));
                        yps_data(counter).yps_signal = mean(cell_curr.signal1);
                        yps_data(counter).cell_pixel_total = length(cell_curr.signal1);
                        yps_data(counter).cell_length_micron = length(cell_curr.signal1)*px_size;
                        yps_data(counter).cell_width_micron = width;
                        yps_data(counter).num_spots = sum(points_in_cell);
                        counter = counter + 1;
                    end
                end
            end
        end
end
% save the dataframe in Excel for sharing the data more broadly.
yps_data_table = struct2table(yps_data);
writetable(yps_data_table,'excel_outputs/grouped_yps_data.xlsx');   
save('mat_outputs/yps_data.mat','yps_data');
warning('on','all');

%% Section 2. Generate Otsu Thresholds
% Now that all the cells are quantified and archived, let's first look at
% each condition (each replicate is pooled) and see whether we can detect
% two populations of low/high fluorescence. In this case, we do detect
% clean bimodal separations. I calculate an Otsu threshold for each that we
% can use to calculate the fraction of cells without plasmids. 
clear; clc; close all;

figsize = [5,5]; % width x height in inches.
save_plot = 1; % 1 to save, 0 to not save.

%-----------------------------------------------------------------------%
load('mat_outputs/yps_data.mat');

% Calculate the bins that encapsulate all the data. We will use the same
% bins for all conditions when plottingso that the histograms are 
% comparable.
[~,edges] = histcounts(log10([yps_data.yps_signal]),100,'Normalization','Probability');
interval = diff(edges);
bin_centers = edges(1:end-1) + interval(1) ./ 2;

f = figure('Units','Inches','Position',[0 0 figsize],'PaperUnits','inches','PaperPosition',[0 0 figsize],'PaperSize',figsize,'CreateFcn','movegui center');

unique_conditions = unique({yps_data(:).condition});
counter=1;

% we want to separate the data by the inducing condition. 
inducing_conditions = {'noninducing','inducing'};
otsu_struct = struct();
otsu_counter = 1;
counter = 1;
for inducing_condition = inducing_conditions
    
    % subset the dataframes according to the current condition. 
    if strcmp(inducing_condition,'noninducing')
        subsetted_conditions = {unique_conditions{contains(unique_conditions,'noninducing')}};
    else
        subsetted_conditions = {unique_conditions{~contains(unique_conditions,'noninducing')}};
    end

    ax = subplot(2,1,counter);

    % Here we seeparate things out by strain for clear, easy archival.
    % wild type first. 
    wt = [yps_data(matches({yps_data(:).condition},subsetted_conditions{2})).yps_signal];
    
    % pull out the wild type counts with the defined edges. 
    [condition_counts,~] = histcounts(log10(wt),edges,'Normalization','Probability');  

    % See the function Hist_Errors at the end of this script. Boot strap
    % the data 100 times using the defined edges to estimate the histogram
    % error
    errors = Hist_Errors(log10(wt),edges);
    
    % plot wild type.
    wt_seb = shadedErrorBar(bin_centers,condition_counts,errors,'LineProps',{'Color',rgb('DodgerBlue'),'LineWidth',1});
    hold on;

    %now subset pcnB
    pcnB = [yps_data(matches({yps_data(:).condition},subsetted_conditions{1})).yps_signal];
    

    [condition_counts,~] = histcounts(log10(pcnB),edges,'Normalization','Probability');
    
    % we want to use the same Otsu threshold per condition, so here measure
    % the otsu threshold using pcnB data since it is a clear bimodal
    % distribution.
    threshold = otsuthresh(condition_counts);
    threshold_idx = round(threshold*length(condition_counts));
    threshold_value = bin_centers(threshold_idx);

    % archive the otsu thresholds for wildtype and pcnB.
    condition = strcat('wt_',inducing_condition);
    otsu_struct(otsu_counter).condition = condition{:};
    otsu_struct(otsu_counter).OtsuThresh = 10^(threshold_value);
    otsu_counter = otsu_counter + 1;
    condition = strcat('pcnB_',inducing_condition);
    otsu_struct(otsu_counter).condition = condition{:};
    otsu_struct(otsu_counter).OtsuThresh = 10^(threshold_value);
    otsu_counter = otsu_counter + 1;
    
    % measure the error of pcnB and plot. 
    errors = Hist_Errors(log10(pcnB),edges);
    pcnb_seb = shadedErrorBar(bin_centers,condition_counts,errors,'LineProps',{'Color',rgb('DarkOrange'),'LineWidth',1});
    xl = xline(threshold_value,'linewidth',1,'Color',rgb('Black'));

    legend([wt_seb.mainLine,pcnb_seb.mainLine,xl],{'Wild Type','\DeltapcnB',['Otsu threshold: ' num2str(10^(threshold_value),4)]},'Box','Off','location','northwest')


    xlabel('log_{10} mean intensity (AU)');
    ylabel('Probability');
    title(inducing_condition,'Interpreter','None');
    set(gca,'Box','Off','TickDir','Out','LineWidth',1,'FontSize',9,'FontName','Arial','YColor','k','YLim',[0,0.25]);
    counter = counter + 1;
end

% Save the data and plot (if desired)
otsu_table = struct2table(otsu_struct);
writetable(otsu_table,'excel_outputs/otsu_thresholds.xlsx');   
save('mat_outputs/otsu_thresholds.mat','otsu_struct');
if save_plot
    print(['plots/histogram_conditions_pooled.pdf'],'-dpdf');
    print(['plots/histogram_conditions_pooled.png'],'-dpng','-r600');
end

%% Section 3. Zoom in on the distributions and compare them. 
% Above I generated the Otsu thersholds. Here, I will zoom in on the
% xlimits [-3,0] to explicitly look at the histograms of wt vs. dpcnB.
% There, we will compare the shape of the distributions with as KS test. I
% will also export the data for these histograms in this zoomed window for
% Jessica to include in the manuscript. 
clear; clc; close all;

figsize = [5,5]; % width x height in inches.
save_plot = 0; % 1 to save, 0 to not save.

%-----------------------------------------------------------------------%
load('mat_outputs/yps_data.mat');

% Calculate the bins that encapsulate all the data. We will use the same
% bins for all conditions when plotting so that the histograms are 
% comparable.
[~,edges] = histcounts(log10([yps_data.yps_signal]),100,'Normalization','Probability');
% based on the plots in Section 2, zoom into histograms edges greater than
% -3.
edges = edges(edges > -3);
interval = diff(edges);
bin_centers = edges(1:end-1) + interval(1) ./ 2;

% archive the bin_centers in a table for export later.
tbl = table(bin_centers','VariableNames',{'x'});

f = figure('Units','Inches','Position',[0 0 figsize],'PaperUnits','inches','PaperPosition',[0 0 figsize],'PaperSize',figsize,'CreateFcn','movegui center');

unique_conditions = unique({yps_data(:).condition});

counter=1;

% As in Section 2, we will loop through condition.
inducing_conditions = {'noninducing','inducing'};
ks_table = table();
for inducing_condition = inducing_conditions

    % subset according to condition.
    if strcmp(inducing_condition,'noninducing')
        subsetted_conditions = {unique_conditions{contains(unique_conditions,'noninducing')}};
    else
        subsetted_conditions = {unique_conditions{~contains(unique_conditions,'noninducing')}};
    end

    % Find the number of cells in each strain of the condition. 
    sizes = [];
    for condition = subsetted_conditions
        data = [yps_data(matches({yps_data(:).condition},condition)).yps_signal];
        sizes = [sizes; length(data)];
    end

    ax = subplot(2,1,counter);

    % pull out the wildtype data.
    wt = [yps_data(matches({yps_data(:).condition},subsetted_conditions{2})).yps_signal];
    % resample the wildtype data without replacement until it meets the
    % smallest N size. 
    wt = datasample(wt,min(sizes),'Replace',false);
    
    % calculate the counts for the predetermined histogram edges.
    [condition_counts,~] = histcounts(log10(wt),edges,'Normalization','Probability');
    
    % archive the wild type counts in the table.
    tbl(:,end+1) = table(condition_counts');

    % calculate the errors of the histogram from bootstrapping. See the
    % Hist_Errors function at the bottom of this script. 
    errors = Hist_Errors(log10(wt),edges);

    % archive the histogram errors into the table.
    tbl(:,end+1) = table(errors');
    wt_seb = shadedErrorBar(bin_centers,condition_counts,errors,'LineProps',{'Color',rgb('DodgerBlue'),'LineWidth',1});
    hold on;

    % now do the same as above to pcnB data. 
    pcnB = [yps_data(matches({yps_data(:).condition},subsetted_conditions{1})).yps_signal];
    pcnB = datasample(pcnB,min(sizes),'Replace',false);
    [condition_counts,~] = histcounts(log10(pcnB),edges,'Normalization','Probability');
    tbl(:,end+1) = table(condition_counts');
    errors = Hist_Errors(log10(pcnB),edges);
    tbl(:,end+1) = table(errors');
    pcnb_seb = shadedErrorBar(bin_centers,condition_counts,errors,'LineProps',{'Color',rgb('DarkOrange'),'LineWidth',1});
    
    legend([wt_seb.mainLine,pcnb_seb.mainLine],{'Wild Type','\DeltapcnB'},'Box','Off','location','northwest')
    
    % compare the wild type and pcnB distributions via a KS test.
    [reject,p,kstatistic] = kstest2(wt,pcnB);
    ks_table(counter,:) = table(inducing_condition,reject,p,kstatistic,'VariableNames',{'condition','reject','p_value','k_statistic'});
    text(-2.75,0.05,['N = ' num2str(length(wt)) 10 'p = ' num2str(p)]);

    xlabel('log_{10} mean intensity (AU)');
    ylabel('Probability');
    title(inducing_condition,'Interpreter','None');
    set(gca,'Box','Off','TickDir','Out','LineWidth',1,'FontSize',9,'FontName','Arial','YColor','k','YLim',[0,0.25]);
    counter = counter + 1;
end
% save the histogram and KS test tables for reference later. 
tbl.Properties.VariableNames = {'x','wt_noninducing','wt_noninducing_sem','pcnb_noninducing','pcnb_noninducing_sem','wt_inducing','wt_inducing_sem','pcnb_inducing','pcnb_inducing_sem'};
writetable(tbl,'excel_outputs/histogram_lineplots.xlsx')
writetable(ks_table,'excel_outputs/ks_test_statistics.xlsx') 

if save_plot
    print(['plots/histograms_with_errorbars.pdf'],'-dpdf');
    print(['plots/histograms_with_errorbars.png'],'-dpng','-r600');
end

%% Section 4. Label the yps data cells that fall below their Otsu thresholds.
% This section will update the source dataframe with a logical true/false
% vector for easy subsetting later. If true, then the given cell is below
% the Otsu threshold and negative for holding plasmid DNA signal. 
clear; clc; close all;

load('mat_outputs/yps_data.mat');
load('mat_outputs/otsu_thresholds.mat')

for idx = 1:length(yps_data)
    condition_curr = yps_data(idx).condition;
    lin_thresh = otsu_struct(matches({otsu_struct(:).condition},condition_curr)).OtsuThresh;

    % create a new field called "below_otsu" that is true or false
    % depending on whether or not the cell's fluorescent signal falls
    % below the Otsu threshold.
    if yps_data(idx).yps_signal < lin_thresh
        yps_data(idx).below_otsu = true;
    else
        yps_data(idx).below_otsu = false;
    end
end

yps_data_table = struct2table(yps_data);
writetable(yps_data_table,'excel_outputs/grouped_yps_data.xlsx');   
save('mat_outputs/yps_data.mat','yps_data');

yps_data = yps_data([yps_data.below_otsu]');
yps_data_table = struct2table(yps_data);
writetable(yps_data_table,'excel_outputs/grouped_yps_data_all_negative.xlsx');   
save('mat_outputs/yps_data_all_negative.mat','yps_data');

%% Section 5. Separate all replicates out and see how the Otsu threshold sits.
% This section splits all the replicates of each condition up and will test
% where the Otsu threshold sits. It  makes a big plot, and saving takes a
% minute. Turn to 1 at your own risk. 
clear; clc; close all;

figsize = [14,10]; % width x height in inches.
save_plot = 0; % 1 to save the plot, 0 to not save the plot. 

%-----------------------------------------------------------------------%
load('mat_outputs/yps_data.mat');
load('mat_outputs/otsu_thresholds.mat')

% look at how the current condition separates with the pre-determined
% histogram edges
[counts,edges] = histcounts(log10([yps_data.yps_signal]),100,'Normalization','Probability');
interval = diff(edges);
bin_centers = edges(1:end-1) + interval(1) ./ 2;

f = figure('Units','Inches','Position',[0 0 figsize],'PaperUnits','inches','PaperPosition',[0 0 figsize],'PaperSize',figsize,'CreateFcn','movegui center');

unique_conditions = unique({yps_data(:).condition});
counter=1;
for condition = unique_conditions
    
    % subset the data to match the current condition
    condition_substruct = yps_data(matches({yps_data(:).condition},condition));
    % obtain the Otsu threshold for the current condition
    lin_thresh = otsu_struct(matches({otsu_struct(:).condition},condition{:})).OtsuThresh;
    % obtain the separate components of the condition for labeling the plot
    % later.
    name_components = split(condition{:},'_');
    for replicate = 1:1:3
        % Grab the data for the current replicate.
        data = [condition_substruct([condition_substruct.replicate] == replicate).yps_signal];
        % Calculate the fraction of cells without plasmids based on the
        % Otsu threshold
        frac_negative = length(data(data < lin_thresh))/length(data);
        
        %plot.
        ax = subplot(length(unique_conditions),3,counter);
        histogram(log10(data),edges,'Normalization','Probability','FaceColor',rgb('GainsBoro'),'EdgeColor','k','LineWidth',0.8);
        hold on;
        xline(log10(lin_thresh),'linewidth',1,'Color',rgb('FireBrick'));
        legend({'Pooled PI stain data',['Otsu threshold: ' num2str(lin_thresh) 10 '% without plasmids: ' num2str(frac_negative)]},'Box','Off','Location','northwest');
        if strcmp(condition,unique_conditions{end})
            xlabel('log_{10} mean intensity (AU)');
        end

        if replicate == 1
            ylabel([name_components{1} '\_' name_components{2} 10 'Frequency'],'FontWeight','Bold');
        end
        
        if strcmp(condition,unique_conditions{1})
            title(['Replicate-' num2str(replicate)])
        end
        
        set(gca,'Box','Off','TickDir','Out','LineWidth',1,'FontSize',9,'FontName','Arial','YColor','k');
        counter = counter + 1;
    end
end

if save_plot
    print(['plots/histogram_conditions_per_replicate.pdf'],'-dpdf');
    print(['plots/histogram_conditions_per_replicate.png'],'-dpng','-r600');
end
%% Section 6. Use the fraction negative from the Otsu threshold to make a bar chart.
% I want to do a direct comparison between conditions to see what the
% fraciton of cells with no fluorescence looks like between conditions.
clear; clc; close all;

figsize = [4,3]; % width x height in inches
save_plot = 0; % 1 to save the plot, 0 to not save the plot. 

%-------------------------------------------------------------------------%
load('mat_outputs/yps_data.mat');
load('mat_outputs/otsu_thresholds.mat')

unique_conditions = unique({yps_data(:).condition});

% Create an empty matrix to populate in the for loop.
frac_negative = [];
condition_counter = 1;
for condition = unique_conditions
    condition = condition{:};

    %split the condition into components
    name_components = split(condition,'_');
    % remake the condition in a way compatible for plotting.
    unique_conditions{condition_counter} = [name_components{1} '\_' name_components{2}];
    % identify the Otsu threshold for the current condition.
    lin_thresh = otsu_struct(matches({otsu_struct(:).condition},condition)).OtsuThresh;
    % subset the data to match the current condition.
    condition_substruct = yps_data(matches({yps_data(:).condition},condition));
    for replicate = 1:1:max([condition_substruct.replicate])
        % pull out the data for the current replicate.
        data = [condition_substruct([condition_substruct.replicate] == replicate).yps_signal];
        % calculate the fraction of cells with no signal.
        frac_negative_temp = length(data(data < lin_thresh))/length(data);
        % add it to the fraction negative matrix.
        frac_negative(replicate,condition_counter) = frac_negative_temp; 
    end
    condition_counter = condition_counter + 1;
end

f = figure('Units','Inches','Position',[0 0 figsize],'PaperUnits','inches','PaperPosition',[0 0 figsize],'PaperSize',figsize,'CreateFcn','movegui center');
% generate a matrix of arbitrary x values.
xdata = repmat([1,2,3,4],3,1);
% make a bar chart with the mean of the fraction negative values
bar([1,2,3,4],mean(frac_negative),'FaceColor',rgb('GainsBoro'));
hold on;
% next a swarm chart of the fraction negative values.
swarmchart(xdata,frac_negative,20,'filled','MarkerFaceColor',rgb('Black'),'MarkerEdgeColor',rgb('Black'),'XJitterWidth',0.25);
xticklabels(unique_conditions);
ylabel('Fraction with no plasmids');

set(gca,'Box','Off','TickDir','Out','LineWidth',1,'FontSize',9,'FontName','Arial','YColor','k');
if save_plot
    print('plots/frac_no_plasmids.pdf','-dpdf');
end
%% Section 7. Calculate the fraction of negative cells for pooled replicates.
% We want to address the big picture: do wild type and pcnB differ? the
% distributions are clear. Here we are extracting values to match the
% distributions in Section 3 for the manuscript. 
clear; clc; close all;

%-------------------------------------------------------------------------%
load('mat_outputs/yps_data.mat');
load('mat_outputs/otsu_thresholds.mat')

unique_conditions = unique({yps_data(:).condition});

frac_negative = struct();
condition_counter = 1;
inducing_conditions = {'noninducing','inducing'};
for inducing_condition = inducing_conditions
    if strcmp(inducing_condition,'noninducing')
        subsetted_conditions = {unique_conditions{contains(unique_conditions,'noninducing')}};
    else
        subsetted_conditions = {unique_conditions{~contains(unique_conditions,'noninducing')}};
    end
    sizes = [];
    for condition = subsetted_conditions
        data = [yps_data(matches({yps_data(:).condition},condition)).yps_signal];
        sizes = [sizes; length(data)];
    end

    for condition = subsetted_conditions
        condition = condition{:};

        % identify the Otsu threshold for the current condition.
        lin_thresh = otsu_struct(matches({otsu_struct(:).condition},condition)).OtsuThresh;
        % subset the data to match the current condition.
        condition_substruct = yps_data(matches({yps_data(:).condition},condition));
        data = [condition_substruct.yps_signal];
        data = datasample(data,min(sizes),'Replace',false);
        % calculate the fraction of cells with no signal.
        frac_negative_temp = length(data(data < lin_thresh))/length(data);
        % add it to the fraction negative structure.
        frac_negative(condition_counter).condition = condition;
        frac_negative(condition_counter).frac_negative = frac_negative_temp;
        frac_negative(condition_counter).N_frac_negative = length(data(data < lin_thresh));
        frac_negative(condition_counter).N_total = length(data);
        condition_counter = condition_counter + 1;
    end
end

tbl = struct2table(frac_negative);
writetable(tbl,'excel_outputs/fracion_negative_ensemble.xlsx')

%% Section 8. Spot count histograms.
% We want to show two things in these data: if cells have spots, how many
% do they have (shown via binomial distributions), and if cells do not have
% spots, then what is the fraction without spots? Here, I construct
% distributions of the number of plasmid foci as well as save the subsetted
% condition datasets.
clear; clc; close all;

figsize = [6,10]; % width x height in inches.
save_plot = 0; % 1 to save the plot, 0 to not save the plot. 

%-----------------------------------------------------------------------%
load('mat_outputs/yps_data.mat');

%filter all data spot counts greater than 0
yps_data = yps_data([yps_data.num_spots]>0);

% look at how the current condition separates with the pre-determined
% histogram edges
[~,edges] = histcounts([yps_data.num_spots],length(unique([yps_data.num_spots])),'Normalization','Probability');
interval = diff(edges);
bin_centers = edges(1:end-1) + interval(1) ./ 2;

f = figure('Units','Inches','Position',[0 0 figsize],'PaperUnits','inches','PaperPosition',[0 0 figsize],'PaperSize',figsize,'CreateFcn','movegui center');

unique_conditions = unique({yps_data(:).condition});
counter=1;
for condition = unique_conditions
    ax = subplot(length(unique_conditions),1,counter);

    % subset the dataframe to the current condition
    condition_substruct = yps_data(matches({yps_data(:).condition},condition));
    % extract only the number of spots
    data = [yps_data(matches({yps_data(:).condition},condition)).num_spots];
    
    histogram(data,edges,'Normalization','Probability','FaceColor',rgb('GainsBoro'),'EdgeColor','k','LineWidth',0.8)

    if counter == 4
        xlabel('Number of plasmid foci (counts)');
    end
    ylabel('Frequency');
    title(condition,'Interpreter','None')
    set(gca,'Box','Off','TickDir','Out','LineWidth',1,'FontSize',9,'FontName','Arial','YColor','k','YLim',[0,0.3]);
    counter = counter + 1;

    % save the subsetted condition dataframe to excel for Jessica's
    % plotting
    condition_tab = struct2table(condition_substruct);
    writetable(condition_tab,['excel_outputs/' condition{:} '_data_without_0_spots.xlsx']);
end
if save_plot
    print('plots/histograms_greater_than_zero.png','-dpng','-r600')
end
%% Section 9. spot histograms per replicate and condition
% Here we will plot all spots per condition and replicate, including counts
% of 0 spots. I will also save the dataframes that include 0-count spots
% for complete data access. This is just to see how all the replicates
% look.
clear; clc; close all;

figsize = [14,10]; % width x height in inches.
save_plot = 0; % 1 to save the plot, 0 to not save the plot. 

%-----------------------------------------------------------------------%
load('mat_outputs/yps_data.mat');

% look at how the current condition separates with the pre-determined
% histogram edges
[~,edges] = histcounts([yps_data.num_spots],length(unique([yps_data.num_spots])),'Normalization','Probability');
interval = diff(edges);
bin_centers = edges(1:end-1) + interval(1) ./ 2;

f = figure('Units','Inches','Position',[0 0 figsize],'PaperUnits','inches','PaperPosition',[0 0 figsize],'PaperSize',figsize,'CreateFcn','movegui center');

unique_conditions = unique({yps_data(:).condition});
counter=1;
for condition = unique_conditions
    %subset the dataframe to the current condition
    condition_substruct = yps_data(matches({yps_data(:).condition},condition));
    %convert the current subsetted dataframe to a table
    condition_tab = struct2table(condition_substruct);
    %save the table to an excel file.
    writetable(condition_tab,['excel_outputs/' condition{:} '_data.xlsx']);
    
    %plot.
    name_components = split(condition{:},'_');
    for replicate = 1:1:3
        ax = subplot(length(unique_conditions),3,counter);
        data = [condition_substruct([condition_substruct.replicate] == replicate).num_spots];

        histogram(data,edges,'Normalization','Probability','FaceColor',rgb('GainsBoro'),'EdgeColor','k','LineWidth',0.8);
        if strcmp(condition,unique_conditions{end})
            xlabel('Number of plasmid foci (counts)');
        end

        if replicate == 1
            ylabel([name_components{1} '\_' name_components{2} 10 'Frequency'],'FontWeight','Bold');
        end
        
        if strcmp(condition,unique_conditions{1})
            title(['Replicate-' num2str(replicate)])
        end
       
        set(gca,'Box','Off','TickDir','Out','LineWidth',1,'FontSize',9,'FontName','Arial','YColor','k');
        counter = counter + 1;
    end
end
if save_plot
    print(['plots/histogram_spots_per_replicate.pdf'],'-dpdf');
    print(['plots/histogram_spots_per_replicate.png'],'-dpng','-r600');
end
%% Section 10. Check the correlation fluorescence vs. spot detection.
% One way to check the accuracy of ThunderSTORM's spot detection is to see
% how well the fraction of cells with 0 spots correlated with the fraction
% of cells with diminished fluorescent signal by Otsu thresholding. I would
% expect it to be around 1:1 as the representative images are fairly clear.

clear; clc; close all;

figsize = [4,3]; % width x height in inches
save_plot  = 1; % 1 to save the plot, 0 to not save the plot. 

%-------------------------------------------------------------------------%
lin_func = @(P,x) P(1)*x+P(2); % linear function, y=mx+b
opts = optimset('Display','off'); % to silence the output of lsqcurvefit
load('mat_outputs/yps_data.mat');
load('mat_outputs/otsu_thresholds.mat')

unique_conditions = unique({yps_data(:).condition});

f = figure('Units','Inches','Position',[0 0 figsize],'PaperUnits','inches','PaperPosition',[0 0 figsize],'PaperSize',figsize,'CreateFcn','movegui center');

condition_counter = 1;
% create empty matrices to populate in the loop below.
frac_negative = []; frac_zero_spots = []; counts_otsu = []; 
counts_otsu_negative = []; counts_spots_negative = [];
counts_spots = []; condition_for_table = {}; replicate_numbers = [];
% loop through each condition.
for condition = unique_conditions
    condition = condition{:};
    name_components = split(condition,'_');
    
    unique_conditions{condition_counter} = [name_components{1} '\_' name_components{2}];
    % pull out the Otsu threshold for the current condition
    lin_thresh = otsu_struct(matches({otsu_struct(:).condition},condition)).OtsuThresh;
    %subset the data to match the current condition.
    condition_substruct = yps_data(matches({yps_data(:).condition},condition));
    
    % loop through each replicate of the current condition.
    for replicate = 1:1:max([condition_substruct.replicate])
        % record the current replicate for saving later
        replicate_numbers = [replicate_numbers; replicate];
        % record the current condition for saving later.
        condition_for_table = [condition_for_table; condition];
        % subset the current replicate's data
        data = [condition_substruct([condition_substruct.replicate] == replicate).yps_signal];
        % find the fraction of negative cells by Otsu threshold
        frac_negative_temp = length(data(data < lin_thresh))/length(data);
        % record the current fraction of negative cells.
        frac_negative = [frac_negative; frac_negative_temp]; 
        % Record the number of cells that are negative. 
        counts_otsu_negative = [counts_otsu_negative; length(data(data < lin_thresh))];
        counts_otsu = [counts_otsu; length(data)];
        % suubset the current replicate's data by number of spots
        data = [condition_substruct([condition_substruct.replicate] == replicate).num_spots];
        % record the fraction of cells with 0 spots
        frac_zero_spots = [frac_zero_spots; length(data(data == 0))/length(data)];
        % record the number of cells with 0 spots. 
        counts_spots = [counts_spots; length(data)];
        counts_spots_negative = [counts_spots_negative; length(data(data == 0))];
    end
    condition_counter = condition_counter + 1;
end

%plotting.
colors_for_plot = [repmat(rgb('DodgerBlue'),6,1);repmat(rgb('DarkOrange'),6,1)]; % specify colors for Replicates 1 and 2.
scatter(frac_negative,frac_zero_spots,50,colors_for_plot,'filled');
hold on

% make blank plots for making a specific legend. 
h(1) = plot(NaN,NaN,'.r','MarkerSize',50,'Color',rgb('DarkOrange'));
h(2) = plot(NaN,NaN,'.r','MarkerSize',50,'Color',rgb('DodgerBlue'));
h(3) = plot(NaN,NaN,':k','LineWidth',1);
h(4) = plot(NaN,NaN,'-k','LineWidth',1);

%generate fits to data.
xfit = linspace(min([frac_negative]),max([frac_negative]),500);
% fit the data with no upper or lower bounds. 
param = lsqcurvefit(lin_func,[1,1],frac_negative,frac_zero_spots,[-inf,-inf],[inf,inf],opts);
% theoretical fit assuming a perfect correlation.
y_theoretical = lin_func([1,0],xfit);
% actual fit with teh data.
y_actual = lin_func(param,xfit);
plot(xfit,y_theoretical,':r','Color',rgb('Black'),'LineWidth',1);
plot(xfit,y_actual,'-r','Color',rgb('Black'),'LineWidth',1)
legend(h,{'Wild type','\Delta\it{pcnB}','Theoretical perfect correlation','Actual fit'},'Box','Off','location','northwest');

xlabel('Fraction negative by Otsu thresholding')
ylabel('Fraction with 0 spots')
set(gca,'Box','Off','TickDir','Out','LineWidth',1,'FontSize',9,'FontName','Arial','YColor','k');
print('plots/fraction_correlation.png','-dpng','-r600')

% save all the recorded features into an excel file.
headers = {'Condition','Replicate','Fraction_by_Otsu','N_otsu_negative','N_otsu','Fraction_by_Counts','N_spots_negative','N_spots'};
SummaryTable = table(condition_for_table,replicate_numbers,frac_negative,counts_otsu_negative,counts_otsu,frac_zero_spots,counts_spots_negative,counts_spots,'VariableNames',headers);
writetable(SummaryTable,['excel_outputs/fraction_negative.xlsx']);

%% Section 11. Plot the distribution of ParB foci in wt and mutant conditions below the Otsu threshold
% The next control to check our accuracy is to plot the number of detected
% spots for all cells below the Otsu threshold. We would expect that the
% predominant histogram result would be a majority of 0 spots. One
% condition, pcnB_inducing, has a significant occurnce of 1-2 spots. This
% makes sense given the correlation above. That says ThunderSTORM
% over-counts the spots a little.
clear; clc; close all;

figsize = [3,5]; % width x height in inches
save_plot  = 1; % 1 to save the plot, 0 to not save the plot. 

%-------------------------------------------------------------------------%
load('mat_outputs/yps_data.mat');
load('mat_outputs/otsu_thresholds.mat');

% filter to include on cells below the Otsu threshold.
yps_data = yps_data([yps_data.below_otsu]);

% look at how the current condition separates with the pre-determined
% histogram edges
[~,edges] = histcounts([yps_data.num_spots],length(unique([yps_data.num_spots])),'Normalization','Probability');
interval = diff(edges);
bin_centers = edges(1:end-1) + interval(1) ./ 2;

unique_conditions = unique({yps_data(:).condition});
counter = 1;
f = figure('Units','Inches','Position',[0 0 figsize],'PaperUnits','inches','PaperPosition',[0 0 figsize],'PaperSize',figsize,'CreateFcn','movegui center');

for condition = unique_conditions
    condition_substruct = yps_data(matches({yps_data(:).condition},condition));
    ax = subplot(length(unique_conditions),1,counter);
    data = [condition_substruct.num_spots];

    histogram(data,edges,'Normalization','Probability','FaceColor',rgb('GainsBoro'),'EdgeColor','k','LineWidth',0.8);
    text(5,0.6,condition,'FontSize',9,'Interpreter','None');
    set(gca,'Box','Off','TickDir','Out','LineWidth',1,'FontSize',9,'FontName','Arial','YColor','k','YLim',[0,1]);
    if counter == 2
        ylabel('Probability')
    end
    if counter == 4
        xlabel('Number of ParB foci (counts)')
    end
    counter = counter + 1;
end
if save_plot
    print(['plots/histogram_spots_from_below_Otsu.pdf'],'-dpdf');
    print(['plots/histogram_spots_from_below_Otsu.png'],'-dpng','-r600');
end
%% Section 12. Violin plots of cell length. 
% Another control plot we need to make is to show how cell length differs
% between cells that fall below/above the Otsu threshold. Wild type here
% apepars to be shorter when they are below the Otsu threshold. 
clear; clc; close all;
figsize = [7,4]; % width x height in inches
save_plot  = 1  ; % 1 to save the plot, 0 to not save the plot. 

%-------------------------------------------------------------------------%
load('mat_outputs/yps_data.mat');
load('mat_outputs/otsu_thresholds.mat');

f = figure('Units','Inches','Position',[0 0 figsize],'PaperUnits','inches','PaperPosition',[0 0 figsize],'PaperSize',figsize,'CreateFcn','movegui center');
unique_conditions = unique({yps_data(:).condition});
counter = 1;
% create empty stacked_data cell aray for populating. 
stacked_data = {};
for condition = unique_conditions
    condition_substruct = yps_data(matches({yps_data(:).condition},condition));
    % Subset for below/above Otsu
    condition_above_Otsu = condition_substruct(~[condition_substruct.below_otsu]);
    condition_below_Otsu = condition_substruct([condition_substruct.below_otsu]);

    % append the stacked data with cell length
    stacked_data{end+1} = [condition_above_Otsu.cell_length_micron];
    stacked_data{end+1} = [condition_below_Otsu.cell_length_micron];

    name_components = split(condition{:},'_');
    unique_conditions{counter} = strcat(name_components{1},'\_',name_components{2});
    counter = counter + 1;
end
% Use padcat to concatenate uneven arrays. 
stacked_data = padcat(stacked_data{:});

colors_for_plot = repmat([rgb('DodgerBlue');rgb('DarkOrange')],4,1); % specify colors for Replicates 1 and 2.

% I edited Hoffman's violin plot from Matlab stack exchange slightly. 
%   - Removed legend inclusion.
%   - Added back in axis ticks.
%   - Disabled median display.
% Source below.
% https://www.mathworks.com/matlabcentral/fileexchange/45134-violin-plot
violin(stacked_data','edgecolor','k','facecolor',colors_for_plot);
xticks([1.5, 3.5, 5.5, 7.5]);
xticklabels(unique_conditions);
ylabel('Cell Length (\mu m)');

%custom legend.
h(1) = bar(NaN,NaN,'FaceColor',rgb('DodgerBlue'),'EdgeColor','none','FaceAlpha',0.5);
h(2) = bar(NaN,NaN,'FaceColor',rgb('DarkOrange'),'EdgeColor','none','FaceAlpha',0.5);
legend(h,{'Above Otsu','Below Otsu'},'Box','Off');
set(gca,'Box','Off','TickDir','Out','LineWidth',1.2,'FontSize',9,'FontName','Arial','YColor','k','YLim',[0,10]);
print(['plots/cell_length_per_GFP_condition.pdf'],'-dpdf');
print(['plots/cell_length_per_GFP_condition.png'],'-dpng','-r600');

%% Section 13. Violin plots of cell width.
% The final control plot we should make would be showing the cell width of
% cells above/below the Otsu threshold. Here is the complete picture. Wild
% type cells below the Otsu threshold are shorter and more narrow. I'm
% inclined to think they are improper segmentations from Oufti. If
% anything, we have far fewer wild type cells that are GFP-free than our
% estimates say. We are over-estimating this value. So the difference
% between wild type and pcnB is potentially greate than we report. Keep in
% mind that wild type is only between 2-7% as of now.
clear; clc; close all;

figsize = [7,4]; % width x height in inches
save_plot  = 1  ; % 1 to save the plot, 0 to not save the plot. 

%-------------------------------------------------------------------------%
load('mat_outputs/yps_data.mat');
load('mat_outputs/otsu_thresholds.mat');

f = figure('Units','Inches','Position',[0 0 figsize],'PaperUnits','inches','PaperPosition',[0 0 figsize],'PaperSize',figsize,'CreateFcn','movegui center');
unique_conditions = unique({yps_data(:).condition});
counter = 1;
stacked_data = {};
for condition = unique_conditions
    condition_substruct = yps_data(matches({yps_data(:).condition},condition));
    condition_above_Otsu = condition_substruct(~[condition_substruct.below_otsu]);
    condition_below_Otsu = condition_substruct([condition_substruct.below_otsu]);
    stacked_data{end+1} = [condition_above_Otsu.cell_width_micron];
    stacked_data{end+1} = [condition_below_Otsu.cell_width_micron];

    name_components = split(condition{:},'_');
    unique_conditions{counter} = strcat(name_components{1},'\_',name_components{2});
    counter = counter + 1;
end
stacked_data = padcat(stacked_data{:});
stacked_shape = size(stacked_data);

colors_for_plot = repmat([rgb('DodgerBlue');rgb('DarkOrange')],4,1); % specify colors for Replicates 1 and 2.
% I edited Hoffman's violin plot from Matlab stack exchange slightly. 
%   - Removed legend inclusion.
%   - Added back in axis ticks.
%   - Disabled median display.
% Source below.
% https://www.mathworks.com/matlabcentral/fileexchange/45134-violin-plot
violin(stacked_data','edgecolor','k','facecolor',colors_for_plot,'medc','');
xticks([1.5, 3.5, 5.5, 7.5]);
xticklabels(unique_conditions);
ylabel('Cell Length (\mu m)');

h(1) = bar(NaN,NaN,'FaceColor',rgb('DodgerBlue'),'EdgeColor','none','FaceAlpha',0.5);
h(2) = bar(NaN,NaN,'FaceColor',rgb('DarkOrange'),'EdgeColor','none','FaceAlpha',0.5);
legend(h,{'Above Otsu','Below Otsu'},'Box','Off');
set(gca,'Box','Off','TickDir','Out','LineWidth',1.2,'FontSize',9,'FontName','Arial','YColor','k');
print(['plots/cell_width_per_GFP_condition.pdf'],'-dpdf');
print(['plots/cell_width_per_GFP_condition.png'],'-dpng','-r600');
%% Functions for this script.

function Errors=Hist_Errors(Data_array,bins_edges)
% Find the error for determining a histogrom by feeding the target
% histogram edges and the source data. Bootstrap the data 100 times and
% make histograms with the input edges. Return the standard deviation of
% the histogram counts in each bin for the error. Note: standard deviation
% of bootstrapped data is the standard error of the mean.
    [bootsam] = bootstrp(100,@(x)x, Data_array);
    store_hist=[];
    
    for i=1:size(bootsam,1)
        [tt1,~]=histcounts(bootsam(i,:),bins_edges,'Normalization','probability');
        store_hist=[store_hist; tt1];
    end
    
    Errors=std(store_hist);
end