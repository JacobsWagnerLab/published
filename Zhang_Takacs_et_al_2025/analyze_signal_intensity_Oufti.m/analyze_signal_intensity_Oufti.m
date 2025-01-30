%% Section 1. Make the main dataframe. 
% Jessica has cells segmented by Oufti, and needs help
% analyzing the number of PI-positive cells out of this dataset. Here, I
% will identify those files, threshold out what is positive and negative,
% then provide analyzed datasets to her. 
% 
% Code constructed by Joshua W. McCausland in the CJW lab, 2024. 
clear; clc; close all;

px_size = 0.065; %microns/px
cell_length_minimum = 5; %minimum cell length for inclusion in the dataset.

% make the subdirectories if they not exist. 
if ~exist('plots', 'dir')
    mkdir('plots')
end
if ~exist('excel_outputs', 'dir')
    mkdir('excel_outputs')
end
if ~exist('mat_outputs', 'dir')
    mkdir('mat_outputs')
end

% Jessica has organized all her files into one directory for me to easily
% access. 
experiment_directory = '/Volumes/Data_02/Jessica Zhang/forJosh_statphasemeshes';

% Identify directories for each biological replicate ("rep1," "rep2")

condition_directories = dir(experiment_directory);
condition_directories = condition_directories(~contains({condition_directories(:).name},'.'));

% This counter is for iteratively appending to the intensity dataframe.
counter = 1;
for condition_index = 1:length(condition_directories)
    current_condition = condition_directories(condition_index).name;
    if ~matches(current_condition,'controls')
        replicate_directories = dir([experiment_directory '/' current_condition '/rep*']); 
        
        for replicate_directory_ind = 1:length(replicate_directories)
            % pull up the current replicate directory and identify all the PI staining
            % Outfi meshes.
            current_replicate_directory = replicate_directories(replicate_directory_ind);
            PI_stain_file_list = dir([current_replicate_directory.folder '/' current_replicate_directory.name '/PI/*.mat']);
            for PI_stain_file_list_ind = 1:length(PI_stain_file_list)
        
                % Identify the current PI stain Oufti mesh. Extract the date from
                % the filename to calculate days for dataframe.
                current_PI_stain_file = PI_stain_file_list(PI_stain_file_list_ind);
                name_components = split(current_PI_stain_file.name,'_');
                date_str = name_components{1};
                date = datetime(date_str,'InputFormat','yyyyMMdd');
        
                if PI_stain_file_list_ind == 1
                    time_0 = date;
                end
                
                %Iterate through every cell and append it to the PI stain
                %dataframe.
                oufti_file_curr = load([current_PI_stain_file.folder '/' current_PI_stain_file.name]);
                for idxa = 1:length(oufti_file_curr.cellList.meshData)
                    region = oufti_file_curr.cellList.meshData{idxa};
                    for idxb = 1:length(region)
                        cell_curr = region{idxb};
                        if length(cell_curr.signal1)*px_size > cell_length_minimum
                            PI_data(counter).date = date;
                            PI_data(counter).condition = current_condition;
                            PI_data(counter).day = days(date-time_0);
                            PI_data(counter).replicate = replicate_directory_ind;
                            PI_data(counter).PI_signal = mean(cell_curr.signal1);
                            PI_data(counter).cell_pixel_total = length(cell_curr.signal1);
                            PI_data(counter).cell_length_micron = length(cell_curr.signal1)*px_size;
                            counter = counter + 1;
                        end
                    end
                end
            end
        end
    else
        controls = dir([experiment_directory '/' current_condition '/*.mat']); 
        for control_ind = 1:length(controls)
            control_name_components = split(controls(control_ind).name,'_');
            control_identity = control_name_components(2);
            oufti_file_curr = load([controls(control_ind).folder '/' controls(control_ind).name]);
            for idxa = 1:length(oufti_file_curr.cellList.meshData)
                region = oufti_file_curr.cellList.meshData{idxa};
                for idxb = 1:length(region)
                    cell_curr = region{idxb};
                    if length(cell_curr.signal1)*px_size > cell_length_minimum
                        PI_data(counter).date = NaN;
                        PI_data(counter).condition = control_identity{:};
                        PI_data(counter).day = NaN;
                        PI_data(counter).replicate = NaN;
                        PI_data(counter).PI_signal = mean(cell_curr.signal1);
                        PI_data(counter).cell_pixel_total = length(cell_curr.signal1);
                        PI_data(counter).cell_length_micron = length(cell_curr.signal1)*px_size;
                        counter = counter + 1;
                    end
                end
            end
        end
    end
end


% save the dataframe in Excel for sharing the data more broadly.
pi_data_table = struct2table(PI_data);
writetable(pi_data_table,'excel_outputs/grouped_pi_data.xlsx');   
save('mat_outputs/pi_data.mat','PI_data');

%% Section 2. Plot all data pooled together.
% This is to view what the data look like before proceeding. 
clear; clc; close all;
% Notes from making pooled plots. 
% B31-MI threshold parameters:
%   - Replicate 2
%   - Compare Day 0 to Day 22
%
% 297 threshold parameters:
%   - Replicate 1
%   - Compare Day 0 to Day 22
%
% K2_pH6 threshold parameters
%   - Use B31-MI's threshold. Barely any positive cells. 
%
% N40 threshold parameters
%   - Replicate 2
%   - Compare Day 2 to Day 11

figsize = [6,6]; % width x height in inches.
replicate = 2; 
strain = 'N40';
save_plots = 1; % 1 if you want to save plots, 0 if not.

%-----------------------------------------------------------------------%
load('mat_outputs/pi_data.mat');
% Calculate the Otsu threshold based on the Day 22 data. It has a clean
% split between negative/positive
if strcmp(strain,'K2_pH6')
    strain = 'B31MI'; %#ok<UNRCH>
end

strain_struct = PI_data(matches({PI_data(:).condition},strain));
f = figure('Units','Inches','Position',[0 0 figsize],'PaperUnits','inches','PaperPosition',[0 0 figsize],'PaperSize',figsize,'CreateFcn','movegui center');

for replicate_idx = 1:1:2
    ax = subplot(2,1,replicate_idx);
    replicate_struct = strain_struct([strain_struct.replicate] == replicate_idx);
    
    if ~strcmp(strain,'N40')
        data = [replicate_struct([replicate_struct.day] == 0 | [replicate_struct.day] == 22).PI_signal]; %#ok<UNRCH>
    else
        data = [replicate_struct([replicate_struct.day] == 2 | [replicate_struct.day] == 11).PI_signal]; %#ok<UNRCH>
    end
    [counts,edges] = histcounts(log10(data),50,'Normalization','Probability');
    interval = diff(edges);
    bin_centers = edges(1:end-1) + interval(1) ./ 2;
    
    threshold = otsuthresh(counts);
    threshold_idx = round(threshold*length(counts));
    threshold_value = bin_centers(threshold_idx);
    
    
    histogram(log10(data),edges,'Normalization','Probability','FaceColor',rgb('GainsBoro'),'EdgeColor','k','LineWidth',0.8);
    hold on;
    xline(threshold_value,'linewidth',1,'Color',rgb('FireBrick'));
    legend({'Pooled PI stain data',['Otsu threshold: ' num2str(10^(threshold_value))]},'Box','Off','Location','northwest');
    xlabel('log_{10} mean intensity (AU)');
    ylabel('Frequency');
    title([strain ' - Replicate ' num2str(replicate_idx)])
    set(gca,'Box','Off','TickDir','Out','LineWidth',1,'FontSize',9,'FontName','Arial','YColor','k');
end
if save_plots
    print(['plots/pooled_histogram_strain-' strain '.pdf'],'-dpdf');
    print(['plots/pooled_histogram_strain-' strain '.png'],'-dpng','-r600');
end

%% Section 3. Compare B31-MI to N40
% The reason why we must pick different Otsu thrsholds is because the
% strains seem to have different intensity distribtuions. This section is
% just to cleanly illustrate just how divergent B31 is from N40. 
%
% This is to view what the data look like before proceeding. 
clear; clc; close all;
% B31-MI threshold parameters:
%   - Replicate 2
%   - Compare Day 0 to Day 22
%
% 297 threshold parameters:
%   - Replicate 1
%   - Compare Day 0 to Day 22
%
% K2_pH6 threshold parameters
%   - Use B31-MI's threshold. Barely any positive cells. 
%
% N40 threshold parameters
%   - Replicate 2
%   - Compare Day 2 to Day 11

figsize = [6,6]; % width x height in inches.
replicate = 2; 
strains = {'B31MI','N40'};
save_plots = 0; % 1 if you want to save plots, 0 if not.

%-----------------------------------------------------------------------%
load('mat_outputs/pi_data.mat');

replicate_struct = PI_data([PI_data.replicate] == replicate);
strain_struct = replicate_struct(matches({replicate_struct(:).condition},strains));
[~,plot_edges] = histcounts(log10([strain_struct.PI_signal]),50,'Normalization','Probability');
interval = diff(plot_edges);
bin_centers = plot_edges(1:end-1) + interval(1) ./ 2;

f = figure('Units','Inches','Position',[0 0 figsize],'PaperUnits','inches','PaperPosition',[0 0 figsize],'PaperSize',figsize,'CreateFcn','movegui center');

for strain_idx = 1:1:2
    ax = subplot(2,1,strain_idx);
    strain_struct = replicate_struct(matches({replicate_struct(:).condition},strains{strain_idx}));
    
    if ~strcmp(strains{strain_idx},'N40')
        data = [strain_struct([strain_struct.day] == 0 | [strain_struct.day] == 22).PI_signal]; %#ok<UNRCH>
    else
        data = [strain_struct([strain_struct.day] == 2 | [strain_struct.day] == 11).PI_signal]; %#ok<UNRCH>
    end
    
    [counts,otsu_edges] = histcounts(log10(data),50,'Normalization','Probability');
    interval = diff(otsu_edges);
    bin_centers = otsu_edges(1:end-1) + interval(1) ./ 2;
    threshold = otsuthresh(counts);
    threshold_idx = round(threshold*length(counts));
    threshold_value = bin_centers(threshold_idx);
    
    
    h1 = histogram(log10(data),plot_edges,'Normalization','Probability','FaceColor',rgb('GainsBoro'),'EdgeColor','k','LineWidth',0.8);
    hold on;
    xline(threshold_value,'linewidth',1,'Color',rgb('FireBrick'));
    legend({'Pooled PI stain data',['Otsu threshold: ' num2str(10^(threshold_value))]},'Box','Off','Location','northwest');
    xlabel('log_{10} mean intensity (AU)');
    ylabel('Frequency');
    title([strains{strain_idx} ' - Replicate ' num2str(replicate)])
    set(gca,'Box','Off','TickDir','Out','LineWidth',1,'FontSize',9,'FontName','Arial','YColor','k');
end
if save_plots
    print(['plots/pooled_histogram_compare-' strains{1} '_' strains{2} '.pdf'],'-dpdf');
    print(['plots/pooled_histogram_compare-' strains{1} '_' strains{2} '.png'],'-dpng','-r600');
end
%% Section 4. Based on the testing above, calculate Otsu for all strains.
% This will iterate through all strains, calculate the Otsu threshold,
% then save the thresholds as Matlab and excel files.
%   - Use the determined replicate and days from section above.
clear; clc; close all;

% load the PI data.
load('mat_outputs/pi_data.mat');

% pull the unique conditions from the dataframe.
unique_conditions = unique({PI_data(:).condition});

% remove the conditions that aren't of the strains of interest. 
unique_conditions(ismember(unique_conditions,{'PIonly','ethanolmesh'})) = [];

for strain_idx = 1:length(unique_conditions)
    strain_curr = unique_conditions(strain_idx);
    strain_curr = strain_curr{:};

    % first filter. If the strain_curr is K2, then change to B31 for
    % threshod determination.
    if strcmp(strain_curr,'K2_pH6')
        strain = 'B31MI'; 
    else
        strain = strain_curr;
    end
    
    % second filter. Only strain 297 uses replicate 1.
    if strcmp(strain,'297')
        replicate = 1;
    else
        replicate = 2;
    end

    % Subset the data to match the current strain.
    strain_struct = PI_data(matches({PI_data(:).condition},strain));
    strain_struct = strain_struct([strain_struct.replicate] == replicate);
    
    % third filter. If the strain is N40, then change which days to query.
    if ~strcmp(strain,'N40')
        data = [strain_struct([strain_struct.day] == 0 | [strain_struct.day] == 22).PI_signal];
    else
        data = [strain_struct([strain_struct.day] == 2 | [strain_struct.day] == 11).PI_signal];
    end

    % generate the Otsu threshold.
    [counts,edges] = histcounts(log10(data),50,'Normalization','Probability');
    interval = diff(edges);
    bin_centers = edges(1:end-1) + interval(1) ./ 2;
    
    threshold = otsuthresh(counts);
    threshold_idx = round(threshold*length(counts));
    threshold_value = bin_centers(threshold_idx);
    lin_thresh = 10^(threshold_value); %convert the log10 bin to linear space.

    % assign the otsu thresholds to a structure for saving later.
    otsu_struct(strain_idx).strain = strain_curr;
    otsu_struct(strain_idx).OtsuThresh = lin_thresh;
end
otsu_table = struct2table(otsu_struct);
writetable(otsu_table,'excel_outputs/otsu_thresholds.xlsx');   
save('mat_outputs/otsu_thresholds.mat','otsu_struct');

%% Section 5. Calculate the fraction of negative PI stain per day
% Construct line plots showing the fraction of cells negative for PI
% staining through time, and save the calculated data to Excel files. This
% applies the Otsu threshold above to determine the fraction.
clear; clc; close all;
strain = 'N40'; %Possible queries: B31MI, K2_pH6, 297, N40
figsize = [3 2.5]; % width x height in inches.
colors_for_plot = [rgb('DodgerBlue');rgb('DarkOrange')]; % specify colors for Replicates 1 and 2.

%-------------------------------------------------------------------------%
load('mat_outputs/pi_data.mat');
load('mat_outputs/otsu_thresholds.mat');

% pull the approriate otsu threshold given the curren strain.
lin_thresh = otsu_struct(matches({otsu_struct(:).strain},strain)).OtsuThresh;

% Subset the primary dataframe according to strian.
strain_struct = PI_data(matches({PI_data(:).condition},strain));

% identify the timecourse days for that strain.
days_under_consideration = rmmissing(unique([strain_struct.day]));

f = figure('Units','Inches','Position',[0 0 figsize],'PaperUnits','inches','PaperPosition',[0 0 figsize],'PaperSize',figsize,'CreateFcn','movegui center');
for rep = 1:1:2
    % subset strain dataset according to current replicate.
    rep_dataset = strain_struct([strain_struct.replicate] == rep);
    frac_negative = []; cell_count = [];
    for day_under_consideration_index = 1:length(days_under_consideration)
        day_under_consideration = days_under_consideration(day_under_consideration_index);
        
        % subset replicate dataset according to the current day.
        data_curr = [rep_dataset([rep_dataset.day] == day_under_consideration).PI_signal];

        %append the fraction negative and cell count vectors with current
        %values.
        frac_negative = [frac_negative; 1-length(data_curr(data_curr > lin_thresh))/length(data_curr)];
        cell_count = [cell_count; length(data_curr)];
    end

    % assign unique identifiers for saving later. 
    if rep == 1
        frac_negative_rep_1 = frac_negative;
        cell_count_rep_1 = cell_count;
    end
    plot(days_under_consideration,frac_negative,'.-','LineWidth',1.2,'MarkerSize',12,'Color',colors_for_plot(rep,:));
    hold on
end
legend({'Replicate 1','Replicate 2'},'Box','Off','Location','southwest');
xlabel('Time (days)','FontSize',9);
ylabel('Fraction PI negative');
set(gca,'Box','Off','TickDir','Out','LineWidth',1,'FontSize',9,'FontName','Arial','YColor','k','YLim',[0,1]);


print(['plots/frac_negative_strain-' strain '.pdf'],'-dpdf');
print(['plots/frac_negative_strain-' strain '.png'],'-dpng','-r600');

% create a data table of the fraction negative and cell counts for saving. 
headers = {'Day','Replicate_1','Cell_Count_1','Replicate_2','Cell_Count_2'};
SummaryTable = table(days_under_consideration',frac_negative_rep_1,cell_count_rep_1,frac_negative,cell_count,'VariableNames',headers);
writetable(SummaryTable,['excel_outputs/fraction_negative_PI_stain_' strain '.xlsx']);

%% Section 6. Histograms of each day for one replicate.
% This is just a broad overview of seeing the distributions for each day of
% a strain/replicate. This gives as sense for how the signal shifts through
% time.
clear; clc; close all;

% editable variables.
strain = 'N40'; %Possible queries: B31MI, K2_pH6, 297, N40
replicate = 2;
save_plots = 0; % 1 if you want to save plots, 0 if not. Saving takes a long time.
figsize = [5 12]; % width x height in inches.

%-------------------------------------------------------------------------%
load('mat_outputs/pi_data.mat');
f = figure('Units','Inches','Position',[0 0 figsize],'PaperUnits','inches','PaperPosition',[0 0 figsize],'PaperSize',figsize,'CreateFcn','movegui center');

%subset the primary dataframe to specific strain.
strain_struct = PI_data(matches({PI_data(:).condition},strain));

% calculate the edges that accomodate all of the data for the strain.
[~,edges] = histcounts(log10([strain_struct.PI_signal]),50);

% Now subset to the specific replicate.
replicate_data = strain_struct([strain_struct.replicate] == replicate);

% identify the days in this timecourse for looping.
days_under_consideration = rmmissing(unique([strain_struct.day]));
for idxa = 1:length(days_under_consideration)
    ax = subplot(length(days_under_consideration),1,idxa);
    
    % Pull the current data in the replicate data. 
    curr_data = [replicate_data([replicate_data.day] == days_under_consideration(idxa)).PI_signal];
    
    % plot that day as a histogram using the determined bin edges above.
    %   - Consistent edges for all histograms allow us to make clean
    %     comparisons between conditions.
    histogram(log10(curr_data),edges,'Normalization','Probability','FaceColor',rgb('GainsBoro'),'EdgeColor','k','LineWidth',0.8);
    
    % figure aesthetics.
    if idxa ~= length(days_under_consideration)
       set(gca,'XColor','None');
    else
        set(gca,'XColor','k');
    end
    xlabel('log_{10} mean intensity (AU)');

    if days_under_consideration(idxa) == 3
        ylabel('Frequency');
    end
    pos_curr = ax.Position;
    ax.Position = [pos_curr(1),pos_curr(2),pos_curr(3),pos_curr(4)*1.25];

    % Determine a consistent position on the axis scale 
    ylimits=get(gca,'ylim');
    ylimit_scale = (ylimits(2)-ylimits(1))./3;

    xlimits = get(gca,'xlim');
    xlimit_scale = (xlimits(2)-xlimits(1))./15;

    % label the current day on the axis. 
    text(xlimit_scale+xlimits(1),ylimits(2)-ylimit_scale,['Day ' num2str(days_under_consideration(idxa))],'FontName','Arial','FontSize',11,'FontWeight','Bold');
end
set(f.Children,'Box','Off','TickDir','Out','LineWidth',1.2,'FontSize',9,'FontName','Arial','YColor','k');

if save_plots
    print(['plots/histograms_by_day_strain-' strain '_rep' num2str(replicate) '.pdf'],'-dpdf');
    print(['plots/histograms_by_day_strain-' strain '_rep' num2str(replicate) '.png'],'-dpng','-r600');
end
%% Section 6. Make a dataframe of the Hoechst data. 
% Jessica prepared separate Oufti mesh files for the Hoechst data. They are
% the exact same cells as the PI staining, just a different fluorescent
% channel.
clear; clc; close all;

px_size = 0.065; %microns/px
cell_length_minimum = 5; %minimum cell length for inclusion in the dataset.

% Jessica has organized all her files into one directory for me to easily
% access. 
experiment_directory = '/Volumes/Data_02/Jessica Zhang/forJosh_statphasemeshes';

% Identify directories for each biological replicate ("rep1," "rep2")

condition_directories = dir(experiment_directory);
condition_directories = condition_directories(~contains({condition_directories(:).name},'.'));

% This counter is for iteratively appending to the intensity dataframe.
counter = 1;
for condition_index = 1:length(condition_directories)
    current_condition = condition_directories(condition_index).name;
    if ~matches(current_condition,'controls')
        replicate_directories = dir([experiment_directory '/' current_condition '/rep*']); 
        
        for replicate_directory_ind = 1:length(replicate_directories)
            % pull up the current replicate directory and identify all the PI staining
            % Outfi meshes.
            current_replicate_directory = replicate_directories(replicate_directory_ind);
            H_stain_file_list = dir([current_replicate_directory.folder '/' current_replicate_directory.name '/Hoechst/*.mat']);
            for H_stain_file_list_ind = 1:length(H_stain_file_list)
        
                % Identify the current PI stain Oufti mesh. Extract the date from
                % the filename to calculate days for dataframe.
                current_H_stain_file = H_stain_file_list(H_stain_file_list_ind);
                name_components = split(current_H_stain_file.name,'_');
                date_str = name_components{1};
                date = datetime(date_str,'InputFormat','yyyyMMdd');
                day = name_components{end};
                day  = str2num(day(2:strfind(name_components{end},'.')));

        
                if H_stain_file_list_ind == 1
                    time_0 = date;
                end
            
                %Iterate through every cell and append it to the PI stain
                %dataframe.
                load([current_H_stain_file.folder '/' current_H_stain_file.name]);
                for idxa = 1:length(cellList.meshData)
                    region = cellList.meshData{idxa};
                    for idxb = 1:length(region)
                        cell_curr = region{idxb};
                        if length(cell_curr.signal1)*px_size > cell_length_minimum
                            H_data(counter).date = date;
                            H_data(counter).condition = current_condition;
                            H_data(counter).day = day;
                            H_data(counter).replicate = replicate_directory_ind;
                            H_data(counter).H_signal = mean(cell_curr.signal1);
                            H_data(counter).cell_pixel_total = length(cell_curr.signal1);
                            H_data(counter).cell_length_micron = length(cell_curr.signal1)*px_size;
                            counter = counter + 1;
                        end
                    end
                end
            end
        end
    end
end
save('mat_outputs/hoescht_data.mat','H_data');


%% Section 7. Compare the Hoechst data between the days.
% This is similar to the PI staining above. First just plot the single cell
% data per day as a histogram. Only considering the B31-MI data. Here, I
% just want to see how the signal shifts through time. 
clear; clc; close all;
figsize = [5 12];

locations = {'northwest','northwest','northwest','northwest','northwest','northwest','northwest','northeast','northeast','northeast'};
f = figure('Units','Inches','Position',[0 0 figsize],'PaperUnits','inches','PaperPosition',[0 0 figsize],'PaperSize',figsize,'CreateFcn','movegui center');
load('mat_outputs/hoescht_data.mat');
H_data = H_data(matches({H_data(:).condition},'B31MI'));
[~,edges] = histcounts(log10([H_data([H_data.H_signal] > 0.02).H_signal]),50);

replicate = 1;

replicate_data = H_data([H_data.replicate] == replicate);
days_under_consideration = rmmissing(unique([replicate_data.day]));

for idxa = 1:length(days_under_consideration)
    ax = subplot(length(days_under_consideration),1,idxa);
    curr_data = [replicate_data([replicate_data.day] == days_under_consideration(idxa)).H_signal];
        
    h1(idxa) = histogram(curr_data(curr_data > 0.02),10.^edges,'Normalization','Probability','FaceColor',rgb('GainsBoro'),'EdgeColor','k','LineWidth',0.8);
    if idxa ~= length(days_under_consideration)
       set(gca,'XColor','None');
    else
        set(gca,'XColor','k');
    end
    xlabel('log_{10} mean intensity (AU)');

    if days_under_consideration(idxa) == 0
        ylabel('Fraction');
    end
    pos_curr = ax.Position;
    ax.Position = [pos_curr(1),pos_curr(2),pos_curr(3),pos_curr(4)*1.25];

    ylimits=get(gca,'ylim');
    ylimit_scale = (ylimits(2)-ylimits(1))./3;
    text(0.02,ylimits(2)-ylimit_scale,['Day ' num2str(days_under_consideration(idxa))],'FontName','Arial','FontSize',11,'FontWeight','Bold');
    
end
edges = h1(1).BinEdges;
interval = diff(edges);
bin_centers = edges(1:end-1) + interval(1) ./ 2;
set(f.Children,'Box','Off','TickDir','Out','LineWidth',1.2,'FontSize',9,'FontName','Arial','YColor','k','XScale','log');
print(['plots/hoescht_histos_with_fit.pdf'],'-dpdf');

%% Section 8. Does the log scale matter?
% Here, we just plot the mean +/- std. deviation across time. Log scale
% mattered for PI data, but as you can see here the difference between
% early and late time points is not that great on a log scale. In fact, it
% looks like a linear decrease on this scale. I think this may be an 
% exponential decay.
clear; clc; close all;
figsize = [5,3];

replicate = 1;
strain = 'B31MI';

f = figure('Units','Inches','Position',[0 0 figsize],'PaperUnits','inches','PaperPosition',[0 0 figsize],'PaperSize',figsize,'CreateFcn','movegui center');
load('mat_outputs/hoescht_data.mat');
replicate_data = H_data([H_data.replicate] == replicate);
days_under_consideration = rmmissing(unique([replicate_data.day]));

stacked_data = {};
mean_vals = [];
std_vals = [];
for idxa = 1:length(days_under_consideration)
    curr_data = [replicate_data([replicate_data.day] == days_under_consideration(idxa)).H_signal].';
    stacked_data{idxa} = curr_data;
    mean_vals = [mean_vals,mean(curr_data)];
    std_vals = [std_vals,std(curr_data)];
end
errorbar(days_under_consideration',mean_vals,std_vals,'o','MarkerSize',5,'Color','k','LineWidth',1)
xlabel('Time (days)')
ylabel('Hoescht intensity (AU)')
set(gca,'Box','Off','TickDir','Out','LineWidth',1.2,'FontSize',9,'FontName','Arial','YColor','k','YScale','log','XLim',[-inf,inf],'YLim',[-inf,inf]);

%% Section 9. Fit the Hoechst data.
% Plotting the Hoechst data on a linear scale revealed it is an exponential
% decay. I fit the data to that function and plot them in this section.
clear; clc; close all;

strain  = 'N40'; %Possible queries: B31MI, K2_pH6, 297, N40
figsize = [5,6]; % width x height in inches.
day_0   = 1; % In case the starting day matters, but it doesn't
save_plot = 0; % 1 if you wish to save the plot, 0 otherwise.

%-------------------------------------------------------------------------%
f = figure('Units','Inches','Position',[0 0 figsize],'PaperUnits','inches','PaperPosition',[0 0 figsize],'PaperSize',figsize,'CreateFcn','movegui center');
load('mat_outputs/hoescht_data.mat');
strain_struct = H_data(matches({H_data(:).condition},strain));
opts = optimset('Display','off'); % This is to silence lsqcurvefit. 
decay_func = @(P,x) P(1)*exp(-x .* P(2))+P(3); %exponential decay function.

params = [];
for replicate = 1:length(unique([strain_struct.replicate]))
    ax = subplot(2,1,replicate);
    replicate_data = strain_struct([strain_struct.replicate] == replicate);
    days_under_consideration = rmmissing(unique([replicate_data.day]));
    
    mean_vals = [];
    std_vals = [];
    for idxa = 1:length(days_under_consideration)
        curr_data = [replicate_data([replicate_data.day] == days_under_consideration(idxa)).H_signal].';
        mean_vals = [mean_vals,mean(curr_data)];
        std_vals = [std_vals,std(curr_data)];
    end
    errorbar(days_under_consideration',mean_vals,std_vals,'o','MarkerSize',5,'Color','k','LineWidth',1);
    hold on
    xlim([-5,25]);
    ylimits = get(gca,'ylim');
    xlimits = get(gca,'xlim');

    % lsqcurve fit for fitting exponential decay function. 
    %   - I define lower and upper bounds with just infinities for an open
    %     variable search.
    %   - Including opts from above to silence the function. 
    %   - lower/upper boundaries necessary to include if we wish to specify
    %     opts.
    starting_assumptions_for_each_term = [1,1,1]; %arbitrary
    lower_bounds = [repmat(-inf,1,3)];
    upper_bounds = [repmat(inf,1,3)];

    param_temp = lsqcurvefit(decay_func,starting_assumptions_for_each_term,days_under_consideration(day_0:end),double(mean_vals(day_0:end)),lower_bounds,upper_bounds,opts);
    params = [params;param_temp]; %in case I wish to review all parameters.

    % Make a line of best fit from the parameters.
    xfit = linspace(days_under_consideration(day_0),days_under_consideration(end),500);
    yfit = decay_func(param_temp,xfit);
    pfit = plot(xfit,yfit,'-r','Color',rgb('FireBrick'),'LineWidth',1);
    xlabel('Time (days in stationary phase)');
    ylabel('Hoescht intensity (AU)');
    title([strain ' - Replicate ' num2str(replicate)]);    
    
    % I only need to specify on one axes what the fit funciton is.
    if replicate == 1
        legend(pfit,['Fit to $y=y_0*e^{-\lambda x} + C$'],'location','northeast','box','off','Interpreter','latex');
    end
    
    % add text to show what the DNA halftime is. 
    % axis limit scales allow us to define general positions that apply to
    % all subsequent graphs.
    ylimit_scale = (ylimits(2)-ylimits(1))./3;
    xlimit_scale = (xlimits(2)-xlimits(1))./4;
    text(xlimits(2)-xlimit_scale,ylimits(2) - ylimit_scale,['\tau_{1/2} = ' num2str(log(2)./param_temp(2),3) ' days']);
    set(gca,'Box','Off','TickDir','Out','LineWidth',1.2,'FontSize',9,'FontName','Arial','YColor','k','XLim',xlim,'YLim',[0,inf]);
end
if save_plot
    print(['plots/hoescht_expfit_' strain '.png'],'-dpng','-r600');
    print(['plots/hoescht_expfit_' strain '.pdf'],'-dpdf');
end
% create a data table of the fit parameters.
headers = {'rep','y_0','lambda','t_c','C'};
SummaryTable = table([1:1:replicate]',params(:,1),params(:,2),log(2)./params(:,2),params(:,3),'VariableNames',headers);
writetable(SummaryTable,['excel_outputs/exp_fit_params_' strain '.xlsx']);

%% Section 10. Export Hoescht fit for all strains/replicates for Jessica.
% similar as above, but no plots. This makes a single excel file with all
% fit data (x/y values for plotting) as well as mean +/- std. dev for
% plotting.
clear; clc; close all;

%-------------------------------------------------------------------------%
load('mat_outputs/hoescht_data.mat');

unique_strains = unique({H_data(:).condition});
opts = optimset('Display','off');

% exponential decay function. y0 * e^(x*lambda) + C
decay_func = @(P,x) P(1)*exp(-x .* P(2))+P(3);
fit_table = table();
fit_table_var_names = {};
hoechst_values = table();
hoescht_values_var_names = {};
params = [];

counter = 1;
stacked_data = {};
for strain = unique_strains
    strain_struct = H_data(matches({H_data(:).condition},strain));
    
    for replicate = 1:length(unique([strain_struct.replicate]))
        replicate_data = strain_struct([strain_struct.replicate] == replicate);
        days_under_consideration = rmmissing(unique([replicate_data.day]));
        stacked_data{end+1} = days_under_consideration';
        mean_vals = [];
        std_vals = [];
        for idxa = 1:length(days_under_consideration)
            curr_data = [replicate_data([replicate_data.day] == days_under_consideration(idxa)).H_signal].';
            mean_vals = [mean_vals,mean(curr_data)];
            std_vals = [std_vals,std(curr_data)];
        end
        stacked_data{end+1} = mean_vals';
        stacked_data{end+1} = std_vals';

        hoescht_values_var_names = [hoescht_values_var_names,[strain{:} '_R' num2str(replicate) '_Days']];
        hoescht_values_var_names = [hoescht_values_var_names,[strain{:} '_R' num2str(replicate) '_mean']];
        hoescht_values_var_names = [hoescht_values_var_names,[strain{:} '_R' num2str(replicate) '_std']];
        % lsqcurve fit for fitting exponential decay function. 
        %   - I define lower and upper bounds with just infinities for an open
        %     variable search.
        %   - Including opts to silence the function. 
        %   - lower/upper boundaries necessary to include if we wish to specify
        %     opts.
        starting_assumptions_for_each_term = [1,1,1]; %arbitrary
        lower_bounds = [repmat(-inf,1,3)];
        upper_bounds = [repmat(inf,1,3)];
        
        param_temp = lsqcurvefit(decay_func,starting_assumptions_for_each_term,days_under_consideration,double(mean_vals),lower_bounds,upper_bounds,opts);
        params = [params;param_temp]; %in case I wish to review all parameters.
    
        % Make a line of best fit from the parameters.
        xfit = linspace(days_under_consideration(1),days_under_consideration(end),500);
        yfit = decay_func(param_temp,xfit);
        fit_table(:,end+1) = table(xfit');
        fit_table_var_names = [fit_table_var_names, [strain{:} '_R' num2str(replicate) '_fit_x']];
        fit_table(:,end+1) = table(yfit');
        fit_table_var_names = [fit_table_var_names, [strain{:} '_R' num2str(replicate) '_fit_y']];
    end
end

stacked_data = padcat(stacked_data{:});
hoechst_values = array2table(stacked_data);
hoechst_values.Properties.VariableNames = hoescht_values_var_names;
fit_table.Properties.VariableNames = fit_table_var_names;

writetable(hoechst_values,['excel_outputs/measured_hoescht_values.xlsx']);
writetable(fit_table,['excel_outputs/hoescht_exp_fits.xlsx']);

%% Section 11. Growth curves.
% I'm wondering whether I can correlate any of the half time data with
% growth rates. I take a very simple approach here. Just fit the linear
% phases of growth curves (on a log scale) to a straight line and
% extrapolate the doubling time from that. I then create a structure
% archiving doubling times and half times for plotting in the next section.
clear; clc; close all;
figsize = [8,2.5]; % width x height in inches.
colors_for_plot = [rgb('DodgerBlue');rgb('DarkOrange')]; % specify colors for Replicates 1 and 2.

gc = table2struct(readtable('cell_density_meaurements.xlsx'));

% fitting terms.
lin_func = @(P,x) P(1)*x + P(2);
opts = optimset('Display','off');
lower_bounds = [-inf, -inf];
upper_bounds = [inf, inf];

f = figure('Units','Inches','Position',[0 0 figsize],'PaperUnits','inches','PaperPosition',[0 0 figsize],'PaperSize',figsize,'CreateFcn','movegui center');
strains = unique({gc(:).strain});

params = [];
counter = 1;
for strain_idx = 1:length(strains)
    strain_curr = strains(strain_idx);
    strain_curr = strain_curr{:};
    exp_fit_struct = table2struct(readtable(['excel_outputs/exp_fit_params_' strain_curr '.xlsx']));
    strain_data = gc(matches({gc(:).strain},strain_curr));
    ax = subplot(1,3,strain_idx);
    hold on
    for replicate_idx = 1:1:2
        replicate_data = strain_data([strain_data.replicate] == replicate_idx);
        plot([replicate_data.day],log10([replicate_data.density]),'.--','Color',colors_for_plot(replicate_idx,:),'MarkerSize',10,'LineWidth',1);
        % B31-MI replicate 2 doesn't have enough data points in the linear
        % phase to fit, so I skip over that.
        if ~(strcmp(strain_curr,'B31MI') && replicate_idx == 2)
            fit_data = [[replicate_data(1:3).day]',log10([replicate_data(1:3).density])'];
            % fit the linear phase using the fitting terms above.
            params_temp = lsqcurvefit(lin_func,[1,1],fit_data(:,1),fit_data(:,2),lower_bounds,upper_bounds,opts);
            params(counter).strain = strain_curr; % current strain.
            params(counter).replicate = replicate_idx; % current replicate.
            params(counter).slope = params_temp(1); % Raw slope on the log scale of the linear phase.
            params(counter).intercept = params_temp(2); % Intercept of the fit
            params(counter).dt = log10(2)./params_temp(1)*24; % doubling time.
            params(counter).tc = exp_fit_struct(replicate_idx).t_c; % measured half time (days).
            counter = counter + 1;

            xfit = linspace(replicate_data(1).day,replicate_data(3).day,500);
            yfit = lin_func(params_temp,xfit);
            plot(xfit,yfit,'-r','LineWidth',1.2,'Color',rgb('FireBrick'))
        end
    end
    title(strain_curr)
    if strain_idx == 2
        xlabel('Culture age in stationary phase (days)')
    end
    if strain_idx == 1
        ylabel('Culture density (log10(counts/ml))')
    end
    set(gca,'Box','Off','TickDir','Out','LineWidth',1.2,'FontSize',9,'FontName','Arial','YColor','k');
end
print('plots/all_growth_curves.pdf','-dpdf')

%% Section 12. Correlate the doubling time to half life.
% I do not see a correlation between measured half life to doubling time. 
close all;
figsize = [5,4];
f = figure('Units','Inches','Position',[0 0 figsize],'PaperUnits','inches','PaperPosition',[0 0 figsize],'PaperSize',figsize,'CreateFcn','movegui center');
plot([params.dt],[params.tc],'.k','MarkerSize',20)
xlabel('doubling time')
ylabel('Hoechst half life')
set(gca,'Box','Off','TickDir','Out','LineWidth',1.2,'FontSize',9,'FontName','Arial','YColor','k');