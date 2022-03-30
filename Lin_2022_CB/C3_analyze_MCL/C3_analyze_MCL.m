
% 2022 Mar 10th: Final version. This code is doing: 
% Analyze data on the mother cell lineage.
% Perform EMD analysis for baseline, and wavelet analysis for the signal

path = {};
path.folder = 'O:\Shares\Data_04\Wei-Hsiang Lin\uF_dataset\exponential\MCL\';
path.save = 'O:\Shares\Data_04\Wei-Hsiang Lin\uF_dataset\exponential\MCL\';
file_name_list = {'M9GlcCA_lac', 'M9GlcCA_ackA', 'M9GlcCA_pta', 'M9GlcCA_adhE', 'M9GlcCA_ldhA'};

% =======================================================================

param = {};
param.Rsc_adjust_factor = 2;                % Adjust for R_ATP ratio for exposure time (405nm = 200ms, 488 = 400ms)
param.sensor_adjust_factor = [6.84 3.33];   % Weighting factors for QUEEN[2m] protein amount 
param.exposure_adjust_factor = 200;         % exposure time (msec)
param.frame_interval = 6;                   % frame interval for fluorescence (min)

param.wave_length_scale = 2;  % Threshold for "long-term" dynamics, in the unit of cell cycles
                              % For exmaple, when using wave_length_scale = 2 (cell cycle),
                              % longer timescale are regarded as "long-term" dynamics
                              
% ============================= Script start ===============================

for DB = 1:length(file_name_list)

% (1) Construct MCL

cd(path.folder);
file_name = strcat('MCL_combined_', file_name_list{DB});
lineage_data = load( strcat(file_name, '.mat') ).lineage_data;

MCL_data = {};
flag = 1;
MCL_minimal_length = 250;

for m = 1:length(lineage_data)
    
    % Load data
    time = lineage_data{m}.MC.Rsc.dataS(:,1); 
    cc_data = lineage_data{m}.MCstat.cc_data; 
    Rsc = lineage_data{m}.MC.Rsc.dataS(:,2); 
    sensor = lineage_data{m}.MC.sensor.dataS(:,2);   
    ccT = lineage_data{m}.MCstat.mean_ccT;
    MC = lineage_data{m}.MC;
    
    if ( length(time) > MCL_minimal_length )
    
        % (A) Analyze R-score data
        [Rsc_data] = time_series_analysis(Rsc, param, ccT);
    
       % (B) Analyze sensor data
        [sensor_data] = time_series_analysis(sensor, param, ccT);
        
        MCL_data{flag}.MC = MC;
        MCL_data{flag}.cc_data = cc_data;
        MCL_data{flag}.ccT = ccT;
        MCL_data{flag}.Rsc_data = Rsc_data;
        MCL_data{flag}.sensor_data = sensor_data;
        MCL_data{flag}.dataset_index = lineage_data{m}.dataset_index;
        MCL_data{flag}.dataset_indnum = lineage_data{m}.dataset_indnum;
        MCL_data{flag}.chamber_index =  lineage_data{m}.chamber_index;
    
        flag = flag + 1;
    
    end
    
end

% ====================================================================

MCL_summary  = {};

if ( ~isempty(MCL_data) )
    
% (2) Averaging of power spectrum

% original data result
PS_collect = {};
    
tag1 = 'Rsc_data'; tag2 = 'PS'; tag3 = 'PS';
PS_collect.Rsc = PS_analysis_summary(MCL_data, tag1, tag2, tag3);

tag1 = 'Rsc_data'; tag2 = 'PS'; tag3 = 'SH_mean';
PS_collect.RscSH = PS_analysis_summary(MCL_data, tag1, tag2, tag3);

tag1 = 'sensor_data'; tag2 = 'PS'; tag3 = 'PS';
PS_collect.sensor = PS_analysis_summary(MCL_data, tag1, tag2, tag3);

tag1 = 'sensor_data'; tag2 = 'PS'; tag3 = 'SH_mean';
PS_collect.sensorSH = PS_analysis_summary(MCL_data, tag1, tag2, tag3);

MCL_summary.PS_collect = PS_collect;

% (3) Averaging of power spectrum

WLPS = {};
tag = 'Rsc_data';
[WLPS.Rsc.PS, WLPS.Rsc.PSSH] = get_waveletPS_summary(MCL_data, tag);
tag = 'sensor_data';
[WLPS.sensor.PS, WLPS.sensor.PSSH] = get_waveletPS_summary(MCL_data, tag);

MCL_summary.waveletPS = WLPS;

% (4) EMD, ccT and GR summary

[Rsc_WL, sensor_WL] = get_EMD_stat(MCL_data);
[wavelet_stat] = get_wavelet_stat(lineage_data, MCL_data);
[ccT_stat, GR_stat] = get_ccT_GR_stat(lineage_data);
    
MCL_summary.EMD.Rsc_WL = Rsc_WL;
MCL_summary.EMD.sensor_WL = sensor_WL;
MCL_summary.wavelet_summary = wavelet_stat;
MCL_summary.ccT_stat = ccT_stat;
MCL_summary.GR_stat = GR_stat;

end



% (End) Save the MCL data

cd(path.save);
save_name = strcat('MCL_PS_', file_name_list{DB});
save( save_name, 'MCL_data' );

cd(path.save);
save_name = strcat('MCL_summary_', file_name_list{DB});
save( save_name, 'MCL_summary' );

end


% =======================================================================
    
function [xdata] = time_series_analysis(x0, param, ccT)
    
    xdata = {};
    
    % (1) Perform EMD analysis to get nature oscillatory mode
    [EMD_output] = EMD_analysis(x0, param.frame_interval);
    
    % (2) Construct baseline by using the trend and slow modes. 
    baseline = x0;
    mode_num = 0;
    IMF = EMD_output.IMF;
    wave_length = EMD_output.wave_length(:,5);  % mean wave length of each mode
    
    for j = 1:length(wave_length)
        if ( wave_length(j) < param.wave_length_scale *ccT )
            % IMF mode with short wave length;            
            mode_num = mode_num + 1;
            baseline = baseline - IMF(:,j);
        end
    end
    
    % Detrended x: 
    xDT = x0 - baseline;
    
    time = param.frame_interval *(1:length(x0));
    
    xdata.x0 = x0;
    xdata.xDT = xDT;
    xdata.time = time;
    xdata.EMD = EMD_output;
    xdata.baseline = baseline; 
    xdata.signal_mode_num = mode_num;
    
    
    %figure; plot(Rsc); hold on; plot(baseline); hold off;
    
    % (3-1) Power spectrum analysis (with statistical shuffle test)    
    [PS] = PS_analysis(xDT, param.frame_interval);
    xdata.PS = PS;
    
    % (3-2) Wavelet analysis (with statistical shuffle test)
    frequency_Hz = ( 60/param.frame_interval )/ 3600;
    xdata.wavelet = wavelet_analysis(xDT, frequency_Hz);
    
    % (3-3) Autocorrelation analysis (with statistical shuffle test)    
    xdata.ACF = ACF_analysis_with_CI(xDT);

end

% =======================================================================
   
