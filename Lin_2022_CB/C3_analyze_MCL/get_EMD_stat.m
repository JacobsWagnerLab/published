
% Statistics for EMD data

function [Rsc_WL, sensor_WL] = get_EMD_stat(MCL_data)

num = length(MCL_data);
max_mode_num = 5;

Rsc_WL = {};
Rsc_WL.mean = NaN(max_mode_num, num);
Rsc_WL.std = NaN(max_mode_num ,num);
Rsc_WL.power = NaN(max_mode_num, num);

sensor_WL = {};
sensor_WL.mean = NaN(max_mode_num, num);
sensor_WL.std = NaN(max_mode_num, num);
sensor_WL.power = NaN(max_mode_num, num);



for m = 1:num
    
    WL_Rsc = MCL_data{m}.Rsc_data.EMD.wave_length;        
    M = min(max_mode_num, size(WL_Rsc,1));
    
    Rsc_WL.mean(1:M,m) = WL_Rsc(1:M,5);
    Rsc_WL.std(1:M,m) = WL_Rsc(1:M,6);
    Rsc_WL.power(1:M,m) =WL_Rsc(1:M,7);
    
    WL_sensor = MCL_data{m}.sensor_data.EMD.wave_length;        
    M = min(max_mode_num, size(WL_sensor,1));
    
    sensor_WL.mean(1:M,m) = WL_sensor(1:M,5);
    sensor_WL.std(1:M,m) = WL_sensor(1:M,6);
    sensor_WL.power(1:M,m) = WL_sensor(1:M,7);
    
end

end

