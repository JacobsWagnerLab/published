
% Batch analysis for several datasets. 
% Extract the basic statistics for all cell cycles.
% Export the data as "ATP_stat".


path.load = '/Volumes/Data_04/Wei-Hsiang Lin/uF_dataset/exponential/CC/';

dataset_name = {'M9GlcCA_lac', 'M9GlcCA_ldhA', 'M9GlcCA_adhE', 'M9GlcCA_pta', 'M9GlcCA_ackA', ...
                'M9Glc_lac', 'M9Xyl_lac', 'M9Gly_lac',};

ATP_stat = {};

for DB = 1:length(dataset_name)

    cd(path.load);
    
    cc_array_temp = load( strcat('CC_', dataset_name{DB}) ).cc_array;
    
    [cc_array] = cell_cycle_basic_stat(cc_array_temp);
    
    ATP_stat.(dataset_name{DB}) = analyze_ATP(cc_array, dataset_name{DB});
    
end

% =====================================================================

function [output] = analyze_ATP(cc_arrayN, nametag)

ATP_calibration_coef = [1.44 -4.68];

cc_arrayN = cell_cycle_filter_tag(cc_arrayN, 'GR_stat', 'ccT', [40 600]);
cc_arrayN = cell_cycle_filter_tag(cc_arrayN, 'GR_stat', 'GR_CoD', [0.8 1]);
cc_arrayN = cell_cycle_filter_tag(cc_arrayN, 'Rsc_stat', 'mean', [1 10]);
cc_arrayN = cell_cycle_filter_tag(cc_arrayN, 'Rsc_stat', 'std', [0 5]);
cc_arrayN = cell_cycle_filter_tag(cc_arrayN, 'Rsc_stat', 'maxdiff', [0 10]);

ATP = {};
ATP.mean = cell_cycle_extract_tag(cc_arrayN, 'Rsc_stat', 'mean');
ATP.std = cell_cycle_extract_tag(cc_arrayN, 'Rsc_stat', 'std');

ATP.CV = ATP.std./ATP.mean;

ATP.mean_abs = ATP_calibration_coef(1)*ATP.mean + ATP_calibration_coef(2);
ATP.std_abs = ATP_calibration_coef(1)*ATP.std;
ATP.CV_abs = ATP.std_abs./ATP.mean_abs;

ATP.maxdiff = cell_cycle_extract_tag(cc_arrayN, 'Rsc_stat', 'maxdiff');
ATP.max= cell_cycle_extract_tag(cc_arrayN, 'Rsc_stat', 'max');
ATP.min = cell_cycle_extract_tag(cc_arrayN, 'Rsc_stat', 'min');

sensor = {};
sensor.mean = cell_cycle_extract_tag(cc_arrayN, 'sensor_stat', 'mean');
sensor.std = cell_cycle_extract_tag(cc_arrayN, 'sensor_stat', 'std');
sensor.maxdiff = cell_cycle_extract_tag(cc_arrayN, 'sensor_stat', 'maxdiff');

GR = {};
GR.GR = cell_cycle_extract_tag(cc_arrayN, 'GR_stat', 'GR');
GR.ccT = cell_cycle_extract_tag(cc_arrayN, 'GR_stat', 'ccT');
GR.midT = cell_cycle_extract_tag(cc_arrayN, 'GR_stat', 'midT');

cc_info = {};
cc_info.FOV = cell_cycle_extract_tag(cc_arrayN, 'cc_info', 'FOV');
cc_info.CB = cell_cycle_extract_tag(cc_arrayN, 'cc_info', 'CB');


% 
output = {};
output.ATP = ATP;
output.sensor = sensor;
output.GR = GR;
output.cc_info = cc_info;
output.nametag = nametag;

end

