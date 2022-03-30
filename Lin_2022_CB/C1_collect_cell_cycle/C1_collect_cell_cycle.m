
% 2022 Mar 10th: Final version. This code is doing: 
% (1) Collect cell cycles for following analysis. 
% (2) Append ATP, cell size data for each cell cycle.
% (3) Append data of mother/daughter cells if the data is available. 
% Cell cycles are filtered for the ones that has proper mother and daughters. 

path.load = 'O:\Shares\Data_04\Wei-Hsiang Lin\uF_dataset\exponential\combined\';
path.save = 'O:\Shares\Data_04\Wei-Hsiang Lin\uF_dataset\exponential\CC\';
dataset_nametag = {'M9GlcCA_lac', 'M9GlcCA_pta', 'M9GlcCA_ackA', 'M9GlcCA_adhE', 'M9GlcCA_ldhA'};

param = {};
param.Rsc_adjust_factor = 2;                % Adjust for R_ATP ratio for exposure time (405nm = 200ms, 488 = 400ms)
param.sensor_adjust_factor = [6.84 3.33];   % Weighting factors for QUEEN[2m] protein amount 
param.exposure_adjust_factor = 200;         % exposure time (msec)
param.time_unit = 1;                        % minute
param.cc_min_length = 20;                   % minimal cell cycle length (min)

% ============================= Script start ===============================

for DB = 1:length(dataset_nametag)

cd(path.load);
load_name = strcat('combined_', dataset_nametag{DB},'.mat');
temp = load(load_name).combined_dataset;
cc_ensemble = temp.cc_ensemble;

% (1) For all tracks in cc_ensemble, collect the cell cycle
%     that has proper division at both ends

[cc_array] = cell_cycle_collection(cc_ensemble, param);

% (2) Perform basic statistics for biosensor level, Rsc and GR.

[cc_array] = cell_cycle_basic_stat(cc_array, param);
[cc_array] = append_Rsc_cc_percentile(cc_array);
[cc_array] = append_Rsc_mother_offspring(cc_array);

% (3) Save files

cd(path.save);
write_name = strcat('CC_', dataset_nametag{DB});
save(write_name, 'cc_array');

end


% ===================================================================

function [cc_array] = cell_cycle_basic_stat(cc_array, param)

Rsc_adjust_factor = param.Rsc_adjust_factor;
sensor_adjust_factor = param.sensor_adjust_factor;
exposure_adjust_factor = param.exposure_adjust_factor;
time_unit = param.time_unit;

% (1) calculate biosensor level and Rsc

for j = 1:length(cc_array)
    
    data = cc_array{j}.data;
    signal = cc_array{j}.signal;
    
    cc_array{j}.Rsc_stat = {};
    cc_array{j}.sensor_stat = {};    
    
    if ( ~isempty(signal) )
    
    Rsc_traj = Rsc_adjust_factor* (signal(:,3)./signal(:,4));
    sensor_traj = sensor_adjust_factor(1)* signal(:,3) + sensor_adjust_factor(2)* signal(:,4);  
    sensor_traj = sensor_traj / exposure_adjust_factor; 
    
    cc_array{j}.Rsc_traj = Rsc_traj;
    cc_array{j}.sensor_traj = sensor_traj;
    
    % Rsc_statistics
    
    Rsc_stat = {};
    
    Rsc_stat.mean = nanmean(Rsc_traj);
    Rsc_stat.std = nanstd(Rsc_traj);
    Rsc_stat.max = max(Rsc_traj);
    Rsc_stat.min = min(Rsc_traj);
    Rsc_stat.maxdiff = max(Rsc_traj) - min(Rsc_traj);
    Rsc_stat.ini = Rsc_traj(1);
    Rsc_stat.end = Rsc_traj(end);
    
    sensor_stat = {};
    
    sensor_stat.mean = nanmean(sensor_traj);
    sensor_stat.std = nanstd(sensor_traj);
    sensor_stat.max = max(sensor_traj);
    sensor_stat.min = min(sensor_traj);
    sensor_stat.maxdiff = max(sensor_traj) - min(sensor_traj);
    
    cc_array{j}.Rsc_stat = Rsc_stat;
    cc_array{j}.sensor_stat = sensor_stat;
    
    end
    
end


% (2) calculate cell size and growth rate

for j = 1:length(cc_array)
    
    data = cc_array{j}.data;
    signal = cc_array{j}.signal;
    
    cc_array{j}.GR_stat = {};
    cc_array{j}.size_stat = {};
    
    if ( ~isempty(signal) )
        
    time = data(:,1);
    size = data(:,3);
    
    GR_stat = {};

    [GR, CoD, fit_size] = GR_slope(time, size, time_unit);
    
    GR_stat.GR = GR;
    GR_stat.GR_CoD = CoD;
    GR_stat.ccT = length(data)-1;
    GR_stat.midT = (time(end)+time(1) )/ 2;
    
    size_stat = {};
    
    size_stat.mean = nanmean(size);
    size_stat.min = min(size);
    size_stat.max = max(size);
    size_stat.ini = size(1);
    size_stat.add1 = size(end)-size(1);
    size_stat.add2 = fit_size(end) - fit_size(1);
    
    cc_array{j}.GR_stat = GR_stat;
    cc_array{j}.size_stat = size_stat;

    end
    
    
end


end

% ===================================================================

function [GR, CoD, fit_size] = GR_slope(time_point, area, time_unit)

log_area = log(area);

% linear fit on log(area)

coef = polyfit(time_point, log_area, 1);
GR = coef(1)/time_unit;
fit_size = exp(coef(1)*time_point + coef(2));

% Calculating coefficient of determination (CoD, R-square)

log_area_mean = mean(log_area);
log_area_est = polyval(coef, time_point); 

SS_tot = sum( (log_area - log_area_mean).^2 );  % Total sum of square
SS_res = sum( (log_area_est - log_area).^2 );   % Residual sum of square

CoD = 1 - (SS_res/SS_tot);

end

% =======================================================================

function [cc_array] = cell_cycle_collection(cc_ensemble, param)

cc_array = {};
flag = 1;

% Script for each chamber

for CB = 1:length(cc_ensemble)

    cc_data = cc_ensemble{CB}.data;
    % [1] frame index
    % [2] object index
    % [3] object area
    % [4] link type
    % [5,6] linked post object(s) (if any)
    % [7] track label
    % [8] mid_y_position
    
    cc_check = cc_ensemble{CB}.track_check;        
    % Column [1]: Check if the beginning type is a proper division
    % 1: proper division
    % 0: has pre-track, but not proper division (brocken track)
    % (-1): no pre-track or track start in frist frame    
    % Column [2]: Check if the ending type is a proper division
    % 1: proper division
    % 0: has pre-track, no proper division
        
    cc_FOV = cc_ensemble{CB}.dataset_indnum;    
    cc_chamber_index = cc_ensemble{CB}.chamber_index;
    cc_note = cc_ensemble{CB}.track_note;    
    
    [div_rec] = get_div_rec(cc_note);
    
    for k = 1:length(cc_check)
        
        if ( cc_check(k,:) == [1 1] )  % "Good" cell cycle
            
            if ( ~isempty(cc_data{k}.data) )  % non-empty
                
                if ( size(cc_data{k}.data,1) > param.cc_min_length) % cc longer than minimal length
                    
                    cc_array{flag} = cc_data{k};
                    cc_array{flag}.div_rec = div_rec(k,:);
                    cc_array{flag}.ID = [cc_FOV CB k];  %FOV, Chamber, cell track ID
                                                            
                    cc_array{flag}.cc_info.FOV = cc_FOV;
                    cc_array{flag}.cc_info.CB = CB;    
                    cc_array{flag}.cc_info.chamber_index = cc_chamber_index;
    
                    flag = flag + 1;
                
                end
            end
            
        end
        
    end
    
    
end


[cc_array] = get_cc_array_lineage_index(cc_array);

end

% =======================================================================

function [div_rec] = get_div_rec(cc_note)

div_rec = zeros(size(cc_note,1), 3);  % [1] mother [2,3] daughter cell index

for k = 1:size(cc_note,1)
    
    mother_ind = cc_note(k,1);  % ind. of mother cell 
    
    if (mother_ind > 0)
        
        div_rec(k,1) = mother_ind;  
        
        if ( div_rec(mother_ind, 2) == 0 )  

            div_rec(mother_ind, 2) = k;  % assign to 1st daughter
        
        elseif( div_rec(mother_ind, 2) > 0 ) 
            
            div_rec(mother_ind, 3) = k;  % assign to 2nd daughter
        end
        
    end
    
end

end

% =======================================================================

function [cc_array] = get_cc_array_lineage_index(cc_array)

% Find cc_array index for mother cell and offsprings

for m = 1:length(cc_array)

ID = cc_array{m}.ID;
div_rec = cc_array{m}.div_rec;
temp_rec = [NaN NaN NaN];  % temporaly index for mother, offspringA, offspringB 

for j = 1:length(cc_array)
    
    dataJ = cc_array{j};
    
    % (1) Test if this is the same FOV and chamber
    ID_test = (ID(1) == dataJ.ID(1)) && (ID(2) == dataJ.ID(2));
    
    if ( ID_test == 1 )
        
        mother_test = ( div_rec(1) == dataJ.ID(3) ) ;
        
        if ( mother_test == 1 )            
            temp_rec(1) = j; 
        end
           
        offspring_testA = ( div_rec(2) == dataJ.ID(3) ); 
        offspring_testB = ( div_rec(3) == dataJ.ID(3) ); 
        
        if ( offspring_testA == 1)
            temp_rec(2) = j;
        end
        if ( offspring_testB == 1)
            temp_rec(3) = j;
        end
                
    end
    
end

cc_array{m}.cc_info.mother_ind_on_cc_array = temp_rec(1);
cc_array{m}.cc_info.offspringA_ind_on_cc_array = temp_rec(2);
cc_array{m}.cc_info.offspringB_ind_on_cc_array = temp_rec(3);

end

end

% =======================================================================