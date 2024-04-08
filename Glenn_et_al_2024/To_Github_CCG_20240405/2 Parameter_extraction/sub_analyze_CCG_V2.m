
% ========================================================================
% Analyze distribution by binning, curve fitting and crossing of curves
% 2024 Feb 
% ========================================================================

function [rec] = sub_analyze_cell_size_distribution_2M(data)

% (1) Generate histogram

pixel2umSQ = 0.0064;  % conversion factor from pixel to um^2
swarmer_ratio = 0.45; % cell area after division

data(:,2) = data(:,2)* pixel2umSQ; 

bin_size = 0.2;
max_value = 4;

[rec_hist, rec_hist_div] = get_histogram(data, bin_size, max_value, swarmer_ratio);


bin_scale = [0:bin_size:max_value]';      % in um^2
bin_mid_scale = bin_scale + (bin_size/2);
finer_scale = 0.01;
finer_rg = (0:finer_scale:max_value)';

C = {};
C{1} = fit_spline(bin_mid_scale, rec_hist(:,1) + rec_hist(:,4), finer_rg);
C{2} = fit_spline(bin_mid_scale, rec_hist(:,2), finer_rg);
C{3} = fit_spline(bin_mid_scale, rec_hist(:,3), finer_rg);
C{4} = fit_spline(bin_mid_scale, rec_hist_div(:,3), finer_rg);

C_array = zeros(size(C{1},1), 4);
for k = 1:4
    C_array(:,k) = C{k}(:,2);
end


% (2) normalize the cell size distribution 

norm_c1 = rec_hist(:,1)/sum(rec_hist(:,1));
norm_c2 = rec_hist(:,2)/sum(rec_hist(:,2));
norm_c3 = rec_hist(:,3)/sum(rec_hist(:,3));
norm_c3_div = rec_hist_div(:,3)/sum(rec_hist_div(:,3));

nC = {};
nC{1} = fit_spline(bin_mid_scale, norm_c1 , finer_rg);
nC{2} = fit_spline(bin_mid_scale, norm_c2 , finer_rg);
nC{3} = fit_spline(bin_mid_scale, norm_c3 , finer_rg);
nC{4} = fit_spline(bin_mid_scale, norm_c3_div, finer_rg);

nC_array = zeros(size(nC{1},1), 4);
for k = 1:4
    nC_array(:,k) = nC{k}(:,2);
end

%

rgA = [1 1.5];
rgB = [1.4 1.8];
rgC = [2 2.5];

size_ini1 = find_crossing(nC_array(:,4), nC_array(:,1), rgA, finer_scale);
size_ini2 = find_crossing(nC_array(:,1), nC_array(:,2), rgB, finer_scale);
size_ini3 = find_crossing(nC_array(:,2), nC_array(:,3), rgC, finer_scale);

%

Area = [size_ini1 size_ini2 size_ini3];  
% swarmer ini, early stalked ini, late staled ini Area

type_count(1) = sum(rec_hist(:,1), 1);   
type_count(2) = sum(rec_hist(:,2), 1); 
type_count(3) = sum(rec_hist(:,3), 1); 

% ========================================================================

rec = {};

rec.hist = rec_hist;
rec.hist_div = rec_hist_div;
rec.Area = Area;
rec.type_count = type_count;

% ========================================================================

figure('position', [1 1 800 250]);

subplot(211);

plot(finer_rg, C_array(:,1), 'b-'); hold on;
plot(bin_mid_scale, (rec_hist(:,1)+rec_hist(:,4)), 'bo'); hold on;

plot(finer_rg, C_array(:,2), 'g-'); hold on;
plot(bin_mid_scale, rec_hist(:,2), 'go'); hold on;

plot(finer_rg, C_array(:,3), 'r-'); hold on;
plot(bin_mid_scale, rec_hist(:,3), 'ro'); hold on;

%plot(C_axis, C_array(:,4), 'r-.'); hold on;
%plot(bin_mid_scale, rec_hist_half(:,3), 'ro'); hold off;

xlim([0 4.5]);

%

subplot(212)

plot(finer_rg, nC_array(:,1), 'b-'); hold on;
plot(bin_mid_scale, norm_c1, 'bo'); hold on;

plot(finer_rg, nC_array(:,2), 'g-'); hold on;
plot(bin_mid_scale, norm_c2, 'go'); hold on;

plot(finer_rg, nC_array(:,3), 'r-'); hold on;
plot(bin_mid_scale, norm_c3, 'ro'); hold on;

plot(finer_rg, nC_array(:,4), 'r-.'); hold on;
plot(bin_mid_scale, norm_c3_div, 'ro'); hold off;

xlim([0 4.5]);
ylim([0 0.5]);

end


% ===================================================================

function [crossing_ind] = find_crossing(x1, x2, rg, finer_scale)

% find crossing position of two curves

rg_ind = rg / finer_scale;

crossing_ind = NaN;

for j = rg_ind(1):rg_ind(2)
    
    temp1 = x1(j) - x2(j);
    temp2 = x1(j+1) - x2(j+1);
    
    if (temp1 > 0) && (temp2 < 0)
        
        crossing_ind = j;
        
        break;
    
    end    

end

crossing_ind = crossing_ind * finer_scale;

end

% ===================================================================

function [output] = fit_spline(x,y,rg)

% fitting data to spline

fitted_data = spline(x,y,rg);

output = [rg fitted_data];

end

% ===================================================================

function [rec_hist, rec_hist_div] = get_histogram(data, bin_size, max_value, swarmer_ratio)

% Generate histogram

max_bin_num = floor(max_value/bin_size) + 1;
max_type = 4;

rec_hist = zeros(max_bin_num, max_type);
rec_hist_div = zeros(max_bin_num, max_type);

% cell size

for j = 1:size(data,1)
    
    ind = floor(data(j,2)/bin_size);
    
    if ( ind <= 0 )
        ind = 1;
    elseif ( ind >= max_bin_num )
        ind = max_bin_num;
    end
    
    if (data(j,1) <= 4)  % checking for bg noise    
        rec_hist(ind, data(j,1)) = rec_hist(ind, data(j,1)) + 1;  
    end
    
end
   
% half-size distribution

for j = 1:size(data,1)
    
    ind = floor(swarmer_ratio * data(j,2)/bin_size);
    
    if ( ind <= 0 )
        ind = 1;
    elseif ( ind >= max_bin_num )
        ind = max_bin_num;
    end
    
    if (data(j,1) <= 4)  % checking for bg noise  
        rec_hist_div(ind, data(j,1)) = rec_hist_div(ind, data(j,1)) + 1;  
    end
    
end

end

% ===================================================================
