
% Manaully load cc_array_type,
% convert into ATP_stat-like data for following correlation analysis

ATP_stat_type = {};

degree = [8 8 8 8 8 8 4 4];

for DB = 1:8
    
    temp = {};
    name_tag = cc_array_type{DB}.dataset_name;
    
    for k = 1:degree(DB)    
        
        if ( ~isempty(cc_array_type{DB}.cc_array{k} ) )
            
            temp{k} = analyze_ATP(cc_array_type{DB}.cc_array{k}, name_tag);    
        
        else
            
            temp{k} = [];

        end
            
    end

    ATP_stat_type{DB}.data = temp;
    ATP_stat_type{DB}.dataset_name = name_tag;
    
end


% =======================================================================

function [output] = analyze_ATP(cc_arrayN, nametag)

cc_arrayN = cell_cycle_filter_tag(cc_arrayN, 'GR_stat', 'ccT', [40 600]);
cc_arrayN = cell_cycle_filter_tag(cc_arrayN, 'GR_stat', 'GR_CoD', [0.8 1]);
cc_arrayN = cell_cycle_filter_tag(cc_arrayN, 'Rsc_stat', 'mean', [1 10]);
cc_arrayN = cell_cycle_filter_tag(cc_arrayN, 'Rsc_stat', 'std', [0 5]);
cc_arrayN = cell_cycle_filter_tag(cc_arrayN, 'Rsc_stat', 'maxdiff', [0 10]);

ATP = {};
ATP.mean = cell_cycle_extract_tag(cc_arrayN, 'Rsc_stat', 'mean');
ATP.std = cell_cycle_extract_tag(cc_arrayN, 'Rsc_stat', 'std');
ATP.maxdiff = cell_cycle_extract_tag(cc_arrayN, 'Rsc_stat', 'maxdiff');
ATP.max= cell_cycle_extract_tag(cc_arrayN, 'Rsc_stat', 'max');
ATP.min = cell_cycle_extract_tag(cc_arrayN, 'Rsc_stat', 'min');
ATP.ini = cell_cycle_extract_tag(cc_arrayN, 'Rsc_stat', 'ini');


ATP_trend = {};
ATP_trend.slopeR_cc = cell_cycle_extract_tag(cc_arrayN, 'Rsc_stat2', 'slope_R_cc_time');
ATP_trend.slopeZ_cc = cell_cycle_extract_tag(cc_arrayN, 'Rsc_stat2', 'slope_Z_cc_time');

sensor = {};
sensor.mean = cell_cycle_extract_tag(cc_arrayN, 'sensor_stat', 'mean');
sensor.std = cell_cycle_extract_tag(cc_arrayN, 'sensor_stat', 'std');
sensor.maxdiff = cell_cycle_extract_tag(cc_arrayN, 'sensor_stat', 'maxdiff');

GR = {};
GR.GR = cell_cycle_extract_tag(cc_arrayN, 'GR_stat', 'GR');
GR.ccT = cell_cycle_extract_tag(cc_arrayN, 'GR_stat', 'ccT');
GR.midT = cell_cycle_extract_tag(cc_arrayN, 'GR_stat', 'midT');

size = {};
size.mean = cell_cycle_extract_tag(cc_arrayN, 'size_stat', 'mean');
size.ini = cell_cycle_extract_tag(cc_arrayN, 'size_stat', 'ini');
size.min = cell_cycle_extract_tag(cc_arrayN, 'size_stat', 'min');
size.add = cell_cycle_extract_tag(cc_arrayN, 'size_stat', 'add1');


cc_info = {};
cc_info.FOV = cell_cycle_extract_tag(cc_arrayN, 'cc_info', 'FOV');
cc_info.CB = cell_cycle_extract_tag(cc_arrayN, 'cc_info', 'CB');

% 

output = {};
output.ATP = ATP;
output.ATP_trend = ATP_trend;
output.sensor = sensor;
output.GR = GR;
output.cc_info = cc_info;
output.nametag = nametag;

end

% =======================================================================

