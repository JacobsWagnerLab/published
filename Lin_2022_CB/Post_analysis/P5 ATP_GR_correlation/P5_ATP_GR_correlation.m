
% Correlation summary plot to excel
% Manaully load "ATP_stat" data structure

dataset_name = {'M9GlcCA_lac', 'M9GlcCA_ldhA', 'M9GlcCA_adhE', 'M9GlcCA_pta', 'M9GlcCA_ackA'};

DB_num = length(dataset_name);
fit_degree = 2;
min_data_point = 10;

bin_data = {};

for DB = 1:DB_num

tag = dataset_name{DB};
xA = ATP_stat.(tag).ATP.mean;
xB = ATP_stat.(tag).ATP.std;
xC = ATP_stat.(tag).ATP_trend.slopeR_cc;
y = ATP_stat.(tag).GR.GR;

% (A) ATP mean
bin_size = 0.25;
bin_endpoints = [2 8];
bin_dataA = x_binning_V5(xA, y, bin_size, bin_endpoints, min_data_point);
fit_dataA = fit_bin_data(bin_dataA.x, bin_dataA.y.mean, fit_degree);

% (B) ATP SD
bin_size = 0.1;
bin_endpoints = [0 3];
bin_dataB = x_binning_V5(xB, y, bin_size, bin_endpoints, min_data_point);
fit_dataB = fit_bin_data(bin_dataB.x, bin_dataB.y.mean, fit_degree);

% (C) ATP slpe
bin_size = 0.5;
bin_endpoints = [-4 4];
bin_dataC = x_binning_V5(xC, y, bin_size, bin_endpoints, min_data_point);
fit_dataC = fit_bin_data(bin_dataC.x, bin_dataC.y.mean, fit_degree);


bin_data{DB}.ATPmean = xA;
bin_data{DB}.ATPstd = xB;
bin_data{DB}.ATPtrend = xC;
bin_data{DB}.GR = y;

bin_data{DB}.A = bin_dataA;
bin_data{DB}.B = bin_dataB;
bin_data{DB}.C = bin_dataC;

bin_data{DB}.Afit = fit_dataA;
bin_data{DB}.Bfit = fit_dataB;
bin_data{DB}.Cfit = fit_dataC;

end

table = {};

LA = length(bin_data{1}.A.x);
LB = length(bin_data{1}.B.x);
LC = length(bin_data{1}.C.x);

table.A.x = bin_data{1}.A.x;
table.B.x = bin_data{1}.B.x;
table.C.x = bin_data{1}.C.x;

table.A.GRmean = NaN(LA, DB_num);
table.A.GRstd = NaN(LA, DB_num);
table.B.GRmean = NaN(LB, DB_num);
table.B.GRstd = NaN(LB, DB_num);
table.C.GRmean = NaN(LC, DB_num);
table.C.GRstd = NaN(LC, DB_num);

for DB = 1:DB_num
    
    table.A.GRmean(:,DB) = bin_data{DB}.A.y.mean;    
    table.A.GRstd(:,DB) = bin_data{DB}.A.y.std; 
    
    table.B.GRmean(:,DB) = bin_data{DB}.B.y.mean;    
    table.B.GRstd(:,DB) = bin_data{DB}.B.y.std; 
    
    table.C.GRmean(:,DB) = bin_data{DB}.C.y.mean;    
    table.C.GRstd(:,DB) = bin_data{DB}.C.y.std; 
    
    
end

% -----------------------------------------------------------

%close all;

colorA = get(gca,'colororder');
colorA(8,:) = [0.2 0.2 0.2];

colorB = [];
colorB(1,:) = colorA(1,:);
colorB(2,:) = [0.2 0.5 0.1];
colorB(3,:) = [0.6 0.9 0.2];
colorB(4,:) = [0.9 0.2 0.55];
colorB(5,:) = [0.3 0.1 0.15];

color1 = colorA;

yRg = [0.001 0.009];

figure('position', [1 1 800 300]);

subplot(131);

for DB = 1:DB_num
    
    x = bin_data{DB}.ATPmean;
    y = bin_data{DB}.GR;
    bin_rec = bin_data{DB}.A;
    fit_rec = bin_data{DB}.Afit;

    xE = fit_rec.x1;
    yE = fit_rec.y1E;
    
    plot(xE, yE, '-', 'color', color1(DB,:)); hold on;
    errorbar(bin_rec.x, bin_rec.y.mean, bin_rec.y.se, 'o', 'color', color1(DB,:)); 
   
    xlim([2.5 7]);
    ylim(yRg);
    
end

subplot(132);

for DB = 1:DB_num
    
    x = bin_data{DB}.ATPstd;
    y = bin_data{DB}.GR;
    bin_rec = bin_data{DB}.B;
    fit_rec = bin_data{DB}.Bfit;

    xE = fit_rec.x1;
    yE = fit_rec.y1E;
    
    plot(xE, yE, '-', 'color', color1(DB,:)); hold on;
    errorbar(bin_rec.x, bin_rec.y.mean, bin_rec.y.se, 'o', 'color', color1(DB,:)); 
   
    xlim([0 1.3]);
    ylim(yRg);
    
end


subplot(133);

for DB = 1:DB_num
    
    x = bin_data{DB}.ATPtrend;
    y = bin_data{DB}.GR;
    bin_rec = bin_data{DB}.C;
    fit_rec = bin_data{DB}.Cfit;

    xE = fit_rec.x1;
    yE = fit_rec.y1E;
    
    plot(xE, yE, '-', 'color', color1(DB,:)); hold on;
    errorbar(bin_rec.x, bin_rec.y.mean, bin_rec.y.se, 'o', 'color', color1(DB,:)); 
   
    xlim([-3.5 3.5]);
    ylim(yRg);
    
end



% =========================================================

function [data] = fit_bin_data(x0, y0, degree)

% polynomial interpolation
%x0 = bin_data{1}.A.x;
%y0 = bin_data{1}.A.y.mean;

temp = [x0 y0];
temp = RemoveNaN(temp,1);

x1 = temp(:,1);
y1 = temp(:,2);

coef = polyfit(x1, y1, degree);
y1E = polyval(coef, x1);

data= {};
data.x0 = x0;
data.y0 = y0;
data.x1 = x1;
data.y1 = y1;
data.degree = degree;
data.coef = coef;
data.y1E = y1E;

end

% =========================================================

