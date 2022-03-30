
% This code generates the data structure for tempogram plot. 
% To plot the data, use another code "plot tempogram".

% ========================================================================
% Load the cell cycle dataset

path = {};
path.load = '/Volumes/Data_04/Wei-Hsiang Lin/uF_dataset/exponential/CC/';

dataset_nametag = {'M9GlcCA_lac', 'M9GlcCA_adhE', 'M9GlcCA_ldhA', 'M9GlcCA_pta' ,'M9GlcCA_ackA',...
                   'M9Glc_lac', 'M9Xyl_lac', 'M9Gly_lac'};

DB_num = length(dataset_nametag);
tempogram_summary = {};

%

for DB = 1:DB_num

cd(path.load);
dataset_name = strcat('CC_', dataset_nametag{DB}, '.mat');   
cc_array = load(dataset_name).cc_array;

cc_arrayN = cc_array;
cc_arrayN = cell_cycle_filter_tag(cc_arrayN, 'GR_stat', 'ccT', [40 600]);
cc_arrayN = cell_cycle_filter_tag(cc_arrayN, 'GR_stat', 'GR_CoD', [0.8 1]);
cc_arrayN = cell_cycle_filter_tag(cc_arrayN, 'Rsc_stat', 'mean', [1 10]);
cc_arrayN = cell_cycle_filter_tag(cc_arrayN, 'Rsc_stat', 'std', [0 5]);
cc_arrayN = cell_cycle_filter_tag(cc_arrayN, 'Rsc_stat', 'maxdiff', [0 10]);

tempogram_summary{DB}.DBname = (dataset_nametag{DB});
tempogram_summary{DB}.data = cc_tempogram(cc_arrayN);

end


% ========================================================================

function [output] = cc_tempogram(cc_array)

% This code use Rsc data on cell cycle (for both original and z-transformed), 
% and interpolate into 0% to 100%

cc_num  = length(cc_array);

traj_matR = NaN(cc_num, 100);
traj_matZ = NaN(cc_num, 100);

mean_vec = NaN(cc_num, 4);    
SD_vec = NaN(cc_num, 4);       
trend_vec = NaN(cc_num, 4);   % interpolated data

for cc = 1:cc_num

    ccdata = cc_array{cc};    
    traj_matR(cc,:) = ccdata.Rsc_stat2.trajR;   
    traj_matZ(cc,:) = ccdata.Rsc_stat2.trajZ; 
    
    mean_vec(cc,1) = ccdata.Rsc_stat.mean;
    SD_vec(cc,1) = ccdata.Rsc_stat.std;

    trend_vec(cc,1) = ccdata.Rsc_stat2.slope_R_cc_time;     
    trend_vec(cc,2) = ccdata.Rsc_stat2.slope_Z_cc_time;
    trend_vec(cc,3) = ccdata.Rsc_stat2.slope_R_abs_time;    
    trend_vec(cc,4) = ccdata.Rsc_stat2.slope_Z_abs_time;
    
end

sorted_matR0 = [trend_vec(:,1) traj_matR];
sorted_matR0 = sortrows(sorted_matR0, 1);
sorted_matR = sorted_matR0(:,2:end);       % Sort all cell cycle with ascending trend

sorted_matR1 = [mean_vec(:,1) traj_matR];
sorted_matR1 = sortrows(sorted_matR1, 1, 'descend');
sorted_matR_mean = sorted_matR1(:,2:end);  % Sort all cell cycle with descending mean

sorted_matR2 = [SD_vec(:,1) traj_matR];
sorted_matR2 = sortrows(sorted_matR2, 1, 'descend');
sorted_matR_SD = sorted_matR2(:,2:end);    % Sort all cell cycle with descending SD

sorted_matZ0 = [trend_vec(:,2) traj_matZ];
sorted_matZ0 = sortrows(sorted_matZ0, 1);
sorted_matZ = sorted_matZ0(:,2:end);       % Sort all cell cycle with ascending trend

% Output data

output= {};
output.sorted_matR = sorted_matR;
output.sorted_matR_mean = sorted_matR_mean;
output.sorted_matR_SD= sorted_matR_SD;
output.sorted_matZ = sorted_matZ;
output.trend_vec = trend_vec;

end

