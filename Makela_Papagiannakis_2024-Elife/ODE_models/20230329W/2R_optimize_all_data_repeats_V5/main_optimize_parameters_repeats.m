
dataset_flag = 3;  % (1) M9GlcCAA, (2) M9Gly, (3) M9Ala

output_array = {};
save_dir = 'F:\Shares\Data_03\Wei-Hsiang Lin\simulation\DNA_limitation_model\20230329W\';

repeatN = 20;

for rp = 1:repeatN

rp 
    
output = {}; 

% (1) Load parameter
[param0] = load_literature_param_V4();
[simudata0] = sub_simulate_single_param(param0{dataset_flag}, dataset_flag);
[output{1}] = sub_calculate_error_V3B(GR_data, AF_SMT_data, simudata0, dataset_flag);
output{1}.param = param0{dataset_flag};


% (2) Multiple_round_optimization

iterate_round = 30;
initial_size = 30000;
initial_range = 0.5;
iterate_size = 3000;
iterate_range = 0.05;

num_vec = iterate_size *ones(iterate_round,1);
num_vec(1) = initial_size;

range_vec = iterate_range* ones(iterate_round,1);
range_vec(1) = initial_range;


for m = 1:length(range_vec)
    
    %m
    
    range_factor = range_vec(m);
    num = num_vec(m);
    param_temp = output{m}.param;
    
    [output{m+1}] = one_round_optimization(range_factor, num, dataset_flag, param_temp, GR_data, AF_SMT_data, simudata0);
    
end

%[param_matrix] = post_analysis(output);
%plot_optimization(output, dataset_flag);

output_array{rp} = output;

end

cd(save_dir);
save('output_array', 'output_array');

%========================================================================

%function plot_optimization(output, dataset_flag)

area_axis = 0.5*(1:24)';  % um^2
xrange = [0 15];

%{
figure;

round_num = size(output,2);

for p = 1:6
    
subplot(2,3,p);

plot(output{1}.rec{p}.area, output{1}.rec{p}.exp, 'b.'); hold on;
plot(output{1}.rec{p}.area, output{1}.rec{p}.simu, 'r*'); hold on;

for m = 2:round_num-1
    
    plot(output{m}.rec{p}.area, output{m}.rec{p}.exp, 'b.'); hold on;
    plot(output{m}.rec{p}.area, output{m}.rec{p}.simu, 'r.'); hold on;

end

plot(output{m}.rec{p}.area, output{m}.rec{p}.exp, 'b.'); hold on;
plot(output{m}.rec{p}.area, output{m}.rec{p}.simu, 'ro'); hold on;

%ylim([0 0.5]);

end
%}

%

%{
figure;

round_num = size(output,2);

for p = 1:3
    
subplot(1,3,p);

ind1 = 2*p-1;
ind2 = 2*p;

plot(output{1}.rec{ind1}.area, output{1}.rec{ind1}.exp, '.', 'color', [0.5 0.5 1]); hold on;
plot(output{1}.rec{ind1}.area, output{1}.rec{ind1}.simu, '.', 'color', [0 0 0.5]); hold on;

plot(output{1}.rec{ind2}.area, output{1}.rec{ind2}.exp, '.', 'color', [1 1 0.5]); hold on;
plot(output{1}.rec{ind2}.area, output{1}.rec{ind2}.simu, '.', 'color', [0.5 0.5 0]); hold on;

for m = 2:round_num-1
    
  %  plot(output{m}.rec{p}.area, output{m}.rec{p}.exp, 'b.'); hold on;
  %  plot(output{m}.rec{p}.area, output{m}.rec{p}.simu, 'r.'); hold on;

end

plot(output{m}.rec{ind1}.area, output{m}.rec{ind1}.exp, '.', 'color', [0.5 0.5 1]); hold on;
plot(output{m}.rec{ind1}.area, output{m}.rec{ind1}.simu, 'o', 'color', [0 0 0.5]); hold on;

plot(output{m}.rec{ind2}.area, output{m}.rec{ind2}.exp, '.', 'color', [1 1 0.5]); hold on;
plot(output{m}.rec{ind2}.area, output{m}.rec{ind2}.simu, 'o', 'color', [0.5 0.5 0]); hold on;

%ylim([0 0.5]);

end
%}

%end

%========================================================================

function [output_opt] = one_round_optimization(range_factor, num, ...
    dataset_flag, param, exp_data, AF_exp_data, simudata)

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
    temp.c  = (1+rtemp(7)) *temp0.c;
    
    %temp.cp = (1+rtemp(8)) *temp0.cp;       
    %temp.xini  = (1+rtemp(9)) *temp0.xini;
    %temp.yini = (1+rtemp(10)) *temp0.yini;    
    %temp.c  = temp0.c;
    
    temp.cp = temp0.cp;    
    temp.xini = temp0.xini;
    temp.yini = temp0.yini;
    temp.Zs = temp0.Zs;
    temp.Vmean = temp0.Vmean;
    
    param_array{j} = temp;
    
end

% Simulate the param before optimization

simudata_0 = sub_simulate_single_param(param, dataset_flag);    
output_0 = sub_calculate_error_V3B(exp_data, AF_exp_data, simudata_0, dataset_flag);
err_ini = output_0.err.max;

% Simulate the multuple parameters in param_array

parfor j = 1:num
    
    simudata_j = sub_simulate_single_param(param_array{j}, dataset_flag);
    output_j = sub_calculate_error_V3B(exp_data, AF_exp_data, simudata_j, dataset_flag);
    
    err_array(j) = output_j.err.max;
    
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
[output_opt] = sub_calculate_error_V3B(exp_data, AF_exp_data, simudata_opt, dataset_flag);
output_opt.param = param_opt;
output_opt.err_array = err_array;


end

%========================================================================

