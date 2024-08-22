
% Extract the output_array results
% Manually load "output_array" for repeated optimization

data_array = output_array;
save_dir = 'F:\Shares\Data_03\Wei-Hsiang Lin\simulation\DNA_limitation_model\20230329W\';

num_repeat = length(data_array);
num_round = length(data_array{1});

opt_array_traj = {};
opt_array_param = {};

for rp = 1:num_repeat

    temp = data_array{rp};
    opt_array_traj{rp} = temp{end}.rec;
    opt_array_param{rp} = temp{end}.param;

end

% (A1) Converting struct to matrix data format

opt_traj = {};  % From 1 to 6. Each one represents a curve in the data.

for type = 1:6
        
    temp_mat1 = [];
    temp_mat2 = [];
    
    for rp = 1:num_repeat           
        
        temp = opt_array_traj{rp}{type};
        [binned_value] = binning(temp);
        
        temp_mat1(:,rp) = binned_value(:,1);
        temp_mat2(:,rp) = binned_value(:,2);
        
    end
    
    opt_traj{type}.area = temp_mat1;
    opt_traj{type}.simu = temp_mat2;
    
end

% (A2) Comparing parameters between repeated optimized results

param_summary = {};
param_summary.ini = get_param(data_array{1}{1}.param);

param_summary.mat = [];

for rp = 1:num_repeat
    
    param_vec = get_param(opt_array_param{rp});
    
    param_summary.mat(:,rp) = param_vec;
    
end

param_summary.mean = mean(param_summary.mat,2);
param_summary.SD = std(param_summary.mat,0,2);


% (A3) Log experimental data points

exp_data = {};

for type = 1:6
    
    temp_area = opt_array_traj{1}{type}.area;
    temp_exp = opt_array_traj{1}{type}.exp;
    
    exp_data{type}.area = temp_area;
    exp_data{type}.exp = temp_exp;
    
end

% Summarize data 

data_summary = {};
data_summary.opt_traj = opt_traj;
data_summary.param_summary = param_summary;
data_summary.exp_data = exp_data;


%======================================================================%

% (B1) Inspect parameters along optimization process

param_summaryB = {};
param_summaryB.ini = get_param(data_array{1}{1}.param);

for rp = 1:num_repeat    
    
    temp_mat = [];
    
    for round = 1:num_round
    
        param_vec = get_param(data_array{rp}{round}.param);
        temp_mat(:,round) = param_vec;
        
    end
    
    param_summaryB.opt{rp} = temp_mat;
    
end

error_stat = {};

error_mat = [];

for rp = 1:num_repeat
    
    for round = 1:num_round
        temp = data_array{rp}{round}.err.max;
        error_mat(rp,round) = temp;      
        
    end
    
end

% (B2) Calculate error-reduce fraction

error_reduction = [];

for rp = 1:num_repeat
    
    err_vec = error_mat(rp,:);
    
    for round = 1:(num_round-1)
        
        diff = err_vec(round+1)-err_vec(round);
        frac = abs(diff)/err_vec(round);
        
        error_reduction(rp,round) = frac;
        
    end
    
end

% Find the curoff round

err_cutoff = 0.01;
err_mean = mean(error_reduction,1);

round_flag = 1;

for round = 1:length(err_mean)
    
    round_flag = round;
    
    if (err_mean(round) < err_cutoff)        
        break;        
    end    
    
end

% Use the parameter of with round number (round_flag+1):
ter_flag = min(round_flag+1, num_round);

opt_array_trajB = {};
opt_array_paramB = {};

for rp = 1:num_repeat

    temp = data_array{rp};
    opt_array_trajB{rp} = temp{ter_flag}.rec;
    opt_array_paramB{rp} = temp{ter_flag}.param;

end

% (B3) Converting struct to matrix data format

opt_trajB = {};  % From 1 to 6. Each one represents a curve in the data.

for type = 1:6
        
    temp_mat1 = [];
    temp_mat2 = [];
    
    for rp = 1:num_repeat           
        
        temp = opt_array_trajB{rp}{type};
        [binned_value] = binning(temp);
        
        temp_mat1(:,rp) = binned_value(:,1);
        temp_mat2(:,rp) = binned_value(:,2);
        
    end
    
    opt_trajB{type}.area = temp_mat1;
    opt_trajB{type}.simu = temp_mat2;
    
    
end

% (B4) Comparing parameters between repeated optimized results

param_summaryB.ini = get_param(data_array{1}{1}.param);

param_summaryB.mat = [];

for rp = 1:num_repeat
    
    param_vec = get_param(opt_array_paramB{rp});
    
    param_summaryB.mat(:,rp) = param_vec;
    
end

param_summaryB.mean = mean(param_summaryB.mat,2);
param_summaryB.SD = std(param_summaryB.mat,0,2);

%=========================================================================
% (C)Summarize data 
 
data_summaryB = {};

data_summaryB.err_cutoff = err_cutoff;
data_summaryB.ter_flag = ter_flag;
data_summaryB.error_mat = error_mat;
data_summaryB.error_reduction = error_reduction;

data_summaryB.opt_trajB = opt_trajB;
data_summaryB.param_summaryB = param_summaryB;
data_summaryB.exp_data = exp_data;

% Repeat simulation using the optimized parameter 
data_summaryC = optimized_simulation(dataset_flag, data_summaryB);

cd(save_dir);
save('data_summaryC', 'data_summaryC');


%=========================================================================

figure;

subplot(121); plot(error_mat', 'o-');
subplot(122); semilogy(error_reduction', 'o-');


% ========================================================== % 

function [binned_value] = binning(temp)

area_range = [1 12];
tick_size = 0.5;
area_ticks = ( area_range(1) :tick_size :area_range(2) )';

cell_area = [];
binned_value = [];

val_temp = [temp.area temp.simu];

for j = 1:size(area_ticks)-1
        
    cell_area(j) = area_ticks(j)+(tick_size/2);
    
    bin_temp = [];    
    
    for k = 1:length(temp.area)
        
        area = val_temp(k,1);
        
        if ( area >= area_ticks(j) ) && ( area < area_ticks(j+1) )
            
            bin_temp = [bin_temp; val_temp(k,:)];
            
        end
        
    end
        
    if (size(bin_temp,1) > 1)
        
        binned_value(j,:) = nanmean(bin_temp,1);
    
    else
        
        binned_value(j,:) = [NaN NaN];

    end
    
    
end

binned_value = remove_NaN_rows(binned_value);

end

%========================================================================

function [output] = remove_NaN_rows(input)

% remove every row that contains NaN

output = [];

for r = 1:size(input,1)
    
    flag = sum( isnan(input(r,:)) );

    if (flag == 0)
        
        output = [output; input(r,:)];

    end
 
end

end

%========================================================================

function [vec] = get_param(input)

    tag = {'r1','r2','K1','K2','d','Zs','c','cp','xini','yini'};
    
    vec = NaN(10,1);
    
    for r = 1:10            
        vec(r) = input.(tag{r});           
    end

end
