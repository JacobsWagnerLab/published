
% Manaully load "ATP_stat" for all cell cycles.
% Create boxplot for the statistics.

dataset_name = {'M9GlcCA_lac',  'M9Glc_lac', 'M9Xyl_lac', 'M9Gly_lac'};
%dataset_name = {'M9GlcCA_lac','M9GlcCA_ldhA', 'M9GlcCA_adhE', 'M9GlcCA_pta', 'M9GlcCA_ackA'};
%dataset_name = {'M9GlcCA_lac','M9GlcCA_lac_lowIPTG'};
      
data = {};

DB_num = length(dataset_name);

for DB = 1:DB_num
    
    DB_name = dataset_name{DB};
    temp = ATP_stat.(DB_name);    
    
    data{DB}.ATPmean = temp.ATP.mean;
    data{DB}.ATPstd = temp.ATP.std;       
    data{DB}.ATPmaxdiff = temp.ATP.maxdiff;
    
    data{DB}.Sensormean = temp.sensor.mean;
    data{DB}.Sensorstd = temp.sensor.std;    
    data{DB}.GR = temp.GR.GR;
    
end

dataU = {};
dataU.ATPmean = [];
dataU.ATPstd = [];
dataU.ATPmaxdiff = [];
dataU.Sensormean = [];
dataU.Sensorstd = [];
dataU.GR = [];

for DB = 1:DB_num
    
    L = length(data{DB}.GR);
    ind_vec = DB*ones(L,1);
    
    temp = [data{DB}.GR ind_vec];
    dataU.GR = [dataU.GR; temp];
    
    temp = [data{DB}.ATPmean ind_vec];
    dataU.ATPmean = [dataU.ATPmean; temp];
    
    temp = [data{DB}.ATPstd ind_vec];
    dataU.ATPstd = [dataU.ATPstd; temp];
    
    temp = [data{DB}.Sensormean ind_vec];
    dataU.Sensormean = [dataU.Sensormean; temp];
    
    temp = [data{DB}.Sensorstd ind_vec];
    dataU.Sensorstd = [dataU.Sensorstd; temp];

end

figure('position', [1 1 500 200]);

subplot(131);
boxplot(dataU.GR(:,1), dataU.GR(:,2) , 'PlotStyle','compact', 'Symbol', '', 'Whisker', 1);
ylim([0 0.012]);

subplot(132);
boxplot(dataU.ATPmean(:,1), dataU.ATPmean(:,2) , 'PlotStyle','compact', 'Symbol', '', 'Whisker', 1);
ylim([3 7]);

subplot(133);
boxplot(dataU.ATPstd(:,1), dataU.ATPstd(:,2) , 'PlotStyle','compact', 'Symbol', '', 'Whisker', 1);
ylim([0 1]);


