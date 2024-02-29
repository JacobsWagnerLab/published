
dataset_flag = 1;  % (1) M9GlcCAA, (2) M9Gly, (3) M9Ala

err_weight = [1 1 1 1 1 1]';
err_weight = err_weight/sum(err_weight);
%err_range = [1 20; 1 10; 1 10];


output=  {};

% (1) Load parameter
[param0] = load_literature_param_V2();
[simudata0] = sub_simulate_single_param(param0{dataset_flag}, dataset_flag);
[output{1}] = sub_calculate_error_V2(GR_data, AF_SMT_data, simudata0, dataset_flag, err_weight);
output{1}.param = param0{dataset_flag};


% (2) Multiple_round_optimization

iterate_round = 2;
initial_size = 100;
initial_range = 0.2;
iterate_size = 50;
iterate_range = 0.1;

num_vec = iterate_size *ones(iterate_round,1);
num_vec(1) = initial_size;

range_vec = iterate_range* ones(iterate_round,1);
range_vec(1) = initial_range;

%num_vec = 500*[5 1 1 1 1 1 1 1 1 1];
%range_vec = 0.1*[3 1 1 1 1 1 1 1 1 1];

for m = 1:length(range_vec)
    
    m
    
    range_factor = range_vec(m);
    num = num_vec(m);
    param_temp = output{m}.param;
    
    [output{m+1}] = one_round_optimization(range_factor, num, ...
        dataset_flag, param_temp, GR_data, AF_SMT_data, simudata0, err_weight);
    
end

[param_matrix] = post_analysis(output);
plot_optimization(output, dataset_flag);

%========================================================================

function plot_optimization(output, dataset_flag)

area_axis = 0.5*(1:24)';  % um^2
xrange = [0 15];

if (dataset_flag == 1)

figure('position', [1 1 1000 300]);

round_num = size(output,2);


subplot(131);

plot(output{1}.nu{1}.axis, output{1}.nu{1}.E, 'bo'); hold on;
plot(output{1}.nu{2}.axis, output{1}.nu{2}.E, 'ro'); hold on;

plot(output{1}.nu{1}.axis, output{1}.nu{1}.S2, 'b--'); hold on;   % WT
plot(output{1}.nu{2}.axis, output{1}.nu{2}.S2, 'r--'); hold on;   % FtsZ

for m = 2:round_num-1
    
    output_m = output{m};
    
    plot(output_m.nu{1}.axis, output_m.nu{1}.S2, 'b:'); hold on;   % WT
    plot(output_m.nu{2}.axis, output_m.nu{2}.S2, 'r:'); hold on;   % FtsZ

end

plot(output{end}.nu{1}.axis, output{end}.nu{1}.S2, 'b-'); hold on;   % WT
plot(output{end}.nu{2}.axis, output{end}.nu{2}.S2, 'r-'); hold on;   % FtsZ

ylim([0 0.5]);

%

subplot(132);
plot(area_axis, output{1}.AF_RNAP.exp{1}, 'bo'); hold on;
plot(area_axis, output{1}.AF_RNAP.exp{2}, 'ro'); hold on;

plot(area_axis, output{1}.AF_RNAP.simu{1}, 'b--'); hold on;   % WT
plot(area_axis, output{1}.AF_RNAP.simu{2}, 'r--'); hold on;   % FtsZ

for m = 2:round_num-1
    
    output_m = output{m};
    
    plot(area_axis, output_m.AF_RNAP.simu{1}, 'b:'); hold on;   % WT
    plot(area_axis, output_m.AF_RNAP.simu{2}, 'r:'); hold on;   % FtsZ

end

plot(area_axis, output{end}.AF_RNAP.simu{1}, 'b-'); hold on;   % WT
plot(area_axis, output{end}.AF_RNAP.simu{2}, 'r-'); hold on;   % FtsZ

xlim(xrange);
ylim([0 1]);

%

subplot(133);
plot(area_axis, output{1}.AF_ribo.exp{1}, 'bo'); hold on;
plot(area_axis, output{1}.AF_ribo.exp{2}, 'ro'); hold on;


plot(area_axis, output{1}.AF_ribo.simu{1}, 'b--'); hold on;   % WT
plot(area_axis, output{1}.AF_ribo.simu{2}, 'r--'); hold on;   % FtsZ
 
for m = 2:round_num
    
    output_m = output{m};
    
    plot(area_axis, output_m.AF_ribo.simu{1}, 'b:'); hold on;   % WT
    plot(area_axis, output_m.AF_ribo.simu{2}, 'r:'); hold on;   % FtsZ

end

plot(area_axis, output{end}.AF_ribo.simu{1}, 'b-'); hold on;   % WT
plot(area_axis, output{end}.AF_ribo.simu{2}, 'r-'); hold on;   % FtsZ
    
xlim(xrange);
ylim([0 1]);



elseif (dataset_flag >1)

figure('position', [1 1 1000 300]);

round_num = size(output,2);

subplot(131);

plot(output{1}.nu{1}.axis, output{1}.nu{1}.E, 'bo'); hold on;
plot(output{1}.nu{2}.axis, output{1}.nu{2}.E, 'ro'); hold on;

plot(output{1}.nu{1}.axis, output{1}.nu{1}.S2, 'b--'); hold on;   % WT
plot(output{1}.nu{2}.axis, output{1}.nu{2}.S2, 'r--'); hold on;   % FtsZ

for m = 2:round_num-1
    
    output_m = output{m};
    
    plot(output_m.nu{1}.axis, output_m.nu{1}.S2, 'b:'); hold on;   % WT
    plot(output_m.nu{2}.axis, output_m.nu{2}.S2, 'r:'); hold on;   % FtsZ

end

plot(output{end}.nu{1}.axis, output{end}.nu{1}.S2, 'b-'); hold on;   % WT
plot(output{end}.nu{2}.axis, output{end}.nu{2}.S2, 'r-'); hold on;   % FtsZ

%xlim([0 0.5]);

%

subplot(132);
plot(area_axis, output{1}.AF_RNAP.exp{1}, 'ko'); hold on;
plot(area_axis, output{1}.AF_RNAP.exp{2}, 'ro'); hold on;

plot(area_axis, output{1}.AF_RNAP.simu{1}, 'k--'); hold on;   % WT
plot(area_axis, output{1}.AF_RNAP.simu{2}, 'r--'); hold on;   % FtsZ

for m = 2:round_num-1
    
    output_m = output{m};
    
    plot(area_axis, output_m.AF_RNAP.simu{1}, 'k:'); hold on;   % WT
    plot(area_axis, output_m.AF_RNAP.simu{2}, 'r:'); hold on;   % FtsZ

end

plot(area_axis, output{end}.AF_RNAP.simu{1}, 'k-'); hold on;   % WT
plot(area_axis, output{end}.AF_RNAP.simu{2}, 'r-'); hold on;   % FtsZ

xlim(xrange);
ylim([0 1]);

%

subplot(133);
plot(area_axis, output{1}.AF_ribo.exp{1}, 'ko'); hold on;
plot(area_axis, output{1}.AF_ribo.exp{2}, 'ro'); hold on;

plot(area_axis, output{1}.AF_ribo.simu{1}, 'k--'); hold on;   % WT
plot(area_axis, output{1}.AF_ribo.simu{2}, 'r--'); hold on;   % FtsZ
 
for m = 2:round_num-1
    
    output_m = output{m};
    
    plot(area_axis, output_m.AF_ribo.simu{1}, 'k:'); hold on;   % WT
    plot(area_axis, output_m.AF_ribo.simu{2}, 'r:'); hold on;   % FtsZ

end

plot(area_axis, output{end}.AF_ribo.simu{1}, 'k-'); hold on;   % WT
plot(area_axis, output{end}.AF_ribo.simu{2}, 'r-'); hold on;   % FtsZ
    
xlim(xrange);
ylim([0 1]);



end

end

%========================================================================

function [param_matrix, err_record] = post_analysis(output);

% Manaully change location to the datasets
round_num = size(output,2);

% (1) Analyze the error minimization

err_record = NaN(round_num, 7);

for m = 1:round_num
    
    temp = output{m}.err.table;

    for j = 1:6
        err_record(m,j) = temp(j);
    end
      
    err_record(m,7) = output{m}.err.sum;
    
end

figure;

for j = 1:7
    
    subplot(3,3,j);
    plot(err_record(:,j), 'o-');

end


% (2) Analyze the error

param_matrix = NaN(10,round_num);

% Rows
% [1] r1, [2] r2, [3] K1, [4] K2, [5] d, [6] Zs, 
% [7] c, [8] cp, [9] xini, [10] yini

% Cols:
% [1,4,5] condition 1 to 3, before optimization
% [2,4,6] condition 1 to 3, after optimization


tag = {'r1','r2','K1','K2','d','Zs','c','cp','xini','yini'};
   
for m = 1:round_num
        
    for r = 1:10 
        param_matrix(r,m) = output{m}.param.(tag{r});
    end
        
end

%xtag = {'ini', 'opt1', 'opt2', 'opt3'};
figure;

for r = 1:10
    
    subplot(2,5,r);
    
    bar(param_matrix(r,:));
    %xticklabels(xtag);
    %ylabel(tag{r});

end

end

%========================================================================

function [output_opt] = one_round_optimization(range_factor, num, ...
    dataset_flag, param, exp_data, AF_exp_data, simudata, err_weight)

%range_factor = 0.1;
%num = 500;

param_array = {};
err_array = zeros(num,1);

for j = 1:num
    
    temp0 = param;
    temp = {};
    
    % generate random variable (uniformly distributed on [param +/- range]
    
    rvec = (2*rand(10,1) - 1);
    rtemp = range_factor* rvec;
    
    temp.r1 = (1+rtemp(1)) *temp0.r1;
    temp.r2 = (1+rtemp(2)) *temp0.r2;
    temp.K1 = (1+rtemp(3)) *temp0.K1;
    temp.K2 = (1+rtemp(4)) *temp0.K2;        
    temp.d  = (1+rtemp(5)) *temp0.d;
    temp.Zs = (1+rtemp(6)) *temp0.Zs;        
    temp.c  = (1+rtemp(7)) *temp0.c;
    temp.cp = (1+rtemp(8)) *temp0.cp;       
    temp.xini  = (1+rtemp(9)) *temp0.xini;
    temp.yini = (1+rtemp(10)) *temp0.yini;
    
    %temp.xini = temp0.xini;
    %temp.yini = temp0.yini;
    
    temp.Vmean = temp0.Vmean;
    
    param_array{j} = temp;
    
end

% Simulate the param before optimization

simudata_0 = sub_simulate_single_param(param, dataset_flag);    
output_0 = sub_calculate_error_V2(exp_data, AF_exp_data, simudata_0, dataset_flag, err_weight);
err_ini = output_0.err.sum;

% Simulate the multuple parameters in param_array

for j = 1:num
    
    simudata_j = sub_simulate_single_param(param_array{j}, dataset_flag);
    output_j = sub_calculate_error_V2(exp_data, AF_exp_data, simudata_j, dataset_flag, err_weight);
    
    err_array(j) = output_j.err.sum;
    
end

% Find the parameter for smallest error
[val,ind] = min(err_array);

if (val < err_ini)
    param_opt = param_array{ind};
else
    param_opt = param;
end
    
% Simulate again using the optimized parameter
[simudata_opt] = sub_simulate_single_param(param_opt, dataset_flag);

% Append original parameter for comparison
simudata_opt.traj.param0 = temp0;

% Plot before and after optimization
[output_opt] = sub_calculate_error_V2(exp_data, AF_exp_data, simudata_opt, dataset_flag, err_weight);
output_opt.param = param_opt;
output_opt.err_array = err_array;


end

%========================================================================

