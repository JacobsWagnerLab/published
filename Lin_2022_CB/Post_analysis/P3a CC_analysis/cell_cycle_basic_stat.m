
function [cc_array] = cell_cycle_basic_stat(cc_array)

% Analyze basic CC statistics

Rsc_adjust_factor = 2;
sensor_adjust_factor = [1 2];
time_unit = 1;

% (1) calculate biosensor level and Rsc

for j = 1:length(cc_array)
    
    data = cc_array{j}.data;
    signal = cc_array{j}.signal;
    
    cc_array{j}.Rsc_stat = {};
    cc_array{j}.sensor_stat = {};
    
    if ( ~isempty(signal) )
    
    Rsc_traj = Rsc_adjust_factor* (signal(:,3)./signal(:,4));
    sensor_traj = sensor_adjust_factor(1)* signal(:,3) + sensor_adjust_factor(2)* signal(:,4);  
    
    cc_array{j}.Rsc_traj = Rsc_traj;
    cc_array{j}.sensor_traj = sensor_traj;
    
    % Rsc_statistics
    
    Rsc_stat = {};
    
    Rsc_stat.mean = nanmean(Rsc_traj);
    Rsc_stat.std = nanstd(Rsc_traj);
    Rsc_stat.max = max(Rsc_traj);
    Rsc_stat.min = min(Rsc_traj);
    Rsc_stat.maxdiff = max(Rsc_traj) - min(Rsc_traj);
    
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

    [GR, CoD] = GR_slope(time, size, time_unit);
    
    GR_stat.GR = GR;
    GR_stat.GR_CoD = CoD;
    GR_stat.ccT = length(data);
    GR_stat.midT = (time(end)+time(1) )/ 2;
    
    size_stat = {};
    
    size_stat.mean = nanmean(size);
    size_stat.min = min(size);
    size_stat.max = max(size);
    size_stat.add1 = size(end)-size(1);
    size_stat.add2 = size(1)* exp(GR_stat.ccT * GR_stat.GR);
    
    cc_array{j}.GR_stat = GR_stat;
    cc_array{j}.size_stat = size_stat;

    end
    
    
end


end


% ===================================================================

function [GR, CoD] = GR_slope(time_point, area, time_unit)

log_area = log(area);

% linear fit on log(area)

coef = polyfit(time_point, log_area, 1);
GR = coef(1)/time_unit;

% Calculating coefficient of determination (CoD, R-square)

log_area_mean = mean(log_area);
log_area_est = polyval(coef, time_point); 

SS_tot = sum( (log_area - log_area_mean).^2 );  % Total sum of square
SS_res = sum( (log_area_est - log_area).^2 );   % Residual sum of square

CoD = 1 - (SS_res/SS_tot);

end

% ===================================================================
