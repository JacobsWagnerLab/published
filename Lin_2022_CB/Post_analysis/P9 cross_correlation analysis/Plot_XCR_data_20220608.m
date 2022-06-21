
% Plot XCR data.
% Manaully load "XCR_data"

% Perform averaging across trajectories

XCR_stat = {};
DB_num = 8;

for DB = 1:DB_num

    data = XCR_data{DB}.data;
    
    Rsc_mat = [];
    sensor_mat = [];
    
    for m = 1:length(data)
        Rsc_mat(m,:) = data{m}.Rsc_GR_xcorr.XCF(:,2);    
        sensor_mat(m,:) = data{m}.sensor_GR_xcorr.XCF(:,2);  
    end
    
    XCR_stat{DB}.Rsc_mat = Rsc_mat;
    XCR_stat{DB}.sensor_mat = sensor_mat;
    
end

for DB = 1:DB_num

    XCR_stat{DB}.Rsc_mean = nanmean(XCR_stat{DB}.Rsc_mat, 1);
    XCR_stat{DB}.Rsc_std = nanstd(XCR_stat{DB}.Rsc_mat, 0, 1);
        
    XCR_stat{DB}.sensor_mean = nanmean(XCR_stat{DB}.sensor_mat, 1);
    XCR_stat{DB}.sensor_std = nanstd(XCR_stat{DB}.sensor_mat, 0, 1);
    
end

% ========================================================================

tau_vec = [-200:1:200]';  % window for calculating XCR

ind_list = [1 6 7 8];
p = length(ind_list);

figure;
yrange1 = [-0.5 0.5];
yrange2 = [-0.5 0.5];
xrange = [-150 150];

for ind = 1:p

DB = ind_list(ind);

subplot(1,p,ind);

%plot(tau_vec, tau_vec*0, 'k-'); hold on;

num = length(XCR_data{DB}.data);

for m = 1:num
    plot(tau_vec, XCR_data{DB}.data{m}.Rsc_GR_xcorr.XCF(:,2), '-', 'color', 0.7*[1 1 1]); hold on;    
end

plot(tau_vec, XCR_stat{DB}.Rsc_mean, 'r-', 'linewidth', 3); hold on;
        
upper_env = XCR_stat{DB}.Rsc_mean + XCR_stat{DB}.Rsc_std;    
lower_env = XCR_stat{DB}.Rsc_mean - XCR_stat{DB}.Rsc_std;
        
plot(tau_vec, upper_env, 'r-.'); hold on;
plot(tau_vec, lower_env, 'r-.'); hold on;
hold off;

ylim(yrange1);    
xlim(xrange);
%grid on;


end


% Comparison (only plot mean)

yrange1 = [-0.3 0.3];
yrange2 = [-0.3 0.3];
xrange = [-150 150];

figure;

for ind = 1:p
   DB = ind_list(ind);
   plot(tau_vec, XCR_stat{DB}.Rsc_mean, '-', 'linewidth', 3); hold on;
end
ylim(yrange1);    
xlim(xrange);
%grid on;


%}

