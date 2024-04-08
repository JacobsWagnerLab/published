%Written by Skye Glenn in the Jacobs-Wagner lab at Stanford University,
%published in 2024
%
%ABOUT
%Section 1: This code first smooths a cell size trajectory then calculates
%instantaneous growth rate (both absolute and normalized for cell size)
%over the complete trajectory. 
%Secton 2: Builds on Section 1. This code includes the input of the time of 
%the appearance of the second MipZ cloud (i.e., the G1/S transition) to 
%calculate G1-specific characteristics.
%Section 3: This code fits an exponential to each cell area trajectory and
%calculates the residuals. It also normalzes the residuals for cell size
%to avoid bias from cells of unusual size or cells with long trajectories. 
%
%INPUTS
%A structure with two fields: time to division (struct.time_to_div) and
%cell area (struct.area). Each column is a cell. I define 2 structures 
%based on cell identity (sw for swarmer progenies, st for stalked 
%progenies).
%The image acquisition interval in minutes (frame_rate).
%The size of window over which to smooth cell size (window).
%If performing DNA replication intation analysis (Section 2), the structure
%must also include a field  for the time of DNA replication (rep_time).
%
%OUTPUTS
%struct.time_norm: normalized cell cycle time between 0 and 1
%struct.avg_time_norm: time_norm shortened by n-1
%struct.area_smooth: smoothed cell area
%struct.diff_growth_rate: instantaneous growth rate
%struct.norm_diff_growth_rate: instantaneous growth rate normalized for
%cell area
%struct.rep_delay: duration of G1 phase
%struct.norm_rep_delay: duration of G1 phase from struct.time_norm
%struct.area_rep: area at the first frame with 2 MipZ clouds (G1/S)
%struct.norm_g1_gr: the average normalized growth rate during G1
%
%DEPENDENCIES
%none
%
%% Define the inputs
frame_rate = 1.5;             %imaging interval in minutes
window = 12;                  %size of the sliding window for smoothing (frames)
sw = swarmer_data_CJW7364;    %abbreviations for the datasets you want
st = stalked_data_CJW7364;    %to work with

%% Section 1
%initialize fields
sw.time_norm = NaN(size(sw.area,1),size(sw.area,2));
sw.avg_time_norm = NaN(size(sw.area,1)-1,size(sw.area,2));
sw.area_smooth = NaN(size(sw.area,1),size(sw.area,2));

%calculate average normalized time and area smoothed over the window
for c = 1:size(sw.area,2)
    sw.time_norm(:,c) = rescale(sw.time_to_div(:,c),0,1);
    for r = 1:sum(~isnan(sw.time_norm(:,c))) - 1
        sw.avg_time_norm(r,c) = (sw.time_norm(r,c) + sw.time_norm(r+1,c))/2;
    end
    sw.area_smooth(1:sum(~isnan(sw.area(:,c))),c) = smooth(sw.area(1:sum(~isnan(sw.area(:,c))),c),window,'moving');
end
%calculate growth rate and growth rate normalized for cell size
sw.diff_growth_rate = diff(sw.area_smooth./frame_rate);
sw.norm_diff_growth_rate = sw.diff_growth_rate./sw.area_smooth(1:(end-1),:);
clear c;

%now the same for stalked cells
st.time_norm = NaN(size(st.area,1),size(st.area,2));
st.avg_time_norm = NaN(size(st.area,1)-1,size(st.area,2));
st.area_smooth = NaN(size(st.area,1),size(st.area,2));

for c = 1:size(st.area,2)
    st.time_norm(:,c) = rescale(st.time_to_div(:,c),0,1);
    for r = 1:sum(~isnan(st.time_norm(:,c))) - 1
        st.avg_time_norm(r,c) = (st.time_norm(r,c) + st.time_norm(r+1,c))/2;
    end
    st.area_smooth(1:sum(~isnan(st.area(:,c))),c) = smooth(st.area(1:sum(~isnan(st.area(:,c))),c),window,'moving');
end
st.diff_growth_rate = diff(st.area_smooth./frame_rate);
st.norm_diff_growth_rate = st.diff_growth_rate./st.area_smooth(1:(end-1),:);
clear c;

%% Section 2: DNA replication
%If also analyzing MipZ segregation, execute the following to find size at 
%replication initiation

sw.rep_delay = abs(sw.time_to_div(1,:)) + sw.rep_time(1,:);

sw.norm_rep_delay = []; sw.area_rep = [];  sw.norm_g1_gr = [];
for c = 1:size(sw.area,2)
    temp_index = find(sw.time_to_div(:,c) == sw.rep_time(1,c));
    sw.norm_rep_delay(1,c) = sw.time_norm(temp_index,c);
    sw.area_rep(1,c) = sw.area(temp_index,c);
    sw.norm_g1_gr(1,c) = mean(sw.norm_diff_growth_rate(1:temp_index,c));
end
clear c; clear temp_index;

%now stalked progenies
st.rep_delay = abs(st.time_to_div(1,:)) + st.rep_time(1,:);

st.norm_rep_delay = []; st.area_rep = [];  st.norm_g1_gr = [];
for c = 1:size(st.area,2)
    temp_index = find(st.time_to_div(:,c) == st.rep_time(1,c));
    st.norm_rep_delay(1,c) = st.time_norm(temp_index,c);
    st.area_rep(1,c) = st.area(temp_index,c);
    st.norm_g1_gr(1,c) = mean(st.norm_diff_growth_rate(1:temp_index,c));
end
clear c; clear temp_index;

%% Section 3
%Execute the below to calculate residuals & residuals normalized for
%cell size, as in Figure 2D

time = sw.time_norm;
char = sw.area_smooth;

residuals = NaN(size(char));
for c=1:size(char,2)
    temp_length = length(char(~isnan(char(:,c))));
    log_char = log(char(1:temp_length,c));
    p = polyfit(time(1:temp_length,c),log_char(:),1);
    y_p = polyval(p,time(1:temp_length,c));
    resid = log_char(:) - y_p;
    residuals(1:temp_length,c) = resid;
    clear temp_length; clear log_char; clear p; clear y_p; clear resid;
end
clear c;

residuals_normalized = NaN(size(char));
for c = 1:size(char,2)
    for r = 1:length(char(~isnan(char(:,c))))
        temp = (residuals(r,c))/char(r,c);
        residuals_normalized(r,c) = temp;
        clear temp;
    end
end
clear c; clear r; clear time; clear char;