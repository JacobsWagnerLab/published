
% This code generate 2-dimensional distribution of ATP and growth rate
% statistics. Manaully load 'ATP_stat' dataset

dataset_name = {'M9GlcCA_lac'};
DB = 1;

ATPmean = ATP_stat.(dataset_name{DB}).ATP.mean;
ATPstd = ATP_stat.(dataset_name{DB}).ATP.std;
GR = ATP_stat.(dataset_name{DB}).GR.GR;

grid = {};
grid.min = [0 2.5];
grid.max = [1.2 6.5];
grid.size = [0.05 0.2];

minimal_sample_size = 5;

[output2] = bin_2D_data_V3(ATPstd, ATPmean, GR, grid, minimal_sample_size);
