
% Load data_summaryC manaully;
save_path = '/Users/wei-hsiang/Desktop/DNA content model/A2_simu_data/0329/';
save_name = 'dataset4B';
dataset = 3;
temp = data_summaryC;

exp_data = temp.exp_data;
opt_traj = temp.opt_trajB;
simudata = temp.simudata;

% Obtain average trajectory over all repeated simulation
[mean_traj, area_axis] = get_mean_trajectory(simudata);

param_summary = temp.param_summaryB;
num_repeat = length(simudata);


% Perform binning
exp_data_binned = {};

w = 0.5;       % area in um^2
xBD = [2 9];   % binninb boundary, area in um^2
    
for m = 1:6

    x = exp_data{m}.area;
    y = exp_data{m}.exp;
    exp_data_binned{m} = binned_average(x,y,w,xBD);

end


% ====================================================================

error_mat = temp.error_mat;
error_reduction = temp.error_reduction;

% ====================================================================

colorB1 = [71 88 248]/255;
colorB2 = [68 52 206]/255;
colorY1 = [251 210 47]/255;
colorY2 = [245 172 80]/255;
colorG1 = [166 168 171]/255;
colorG2 = [48 48 48]/255;

alpha_val = [0.15 0.4 0.4];


% Compare experiments and models

h3B = figure('position',[1 1 800 200]);

for type = 1:3

    subplot(1,3,type);
    
    ind1 = 2*type - 1;
    ind2 = 2*type;
    
    % (1) Plot experimental results
    
    data_area2 = exp_data{ind2}.area;
    data_exp2 = exp_data{ind2}.exp;    
    scatter(data_area2, data_exp2, 6, colorY1, 'filled', 'MarkerFaceAlpha', alpha_val(type), 'MarkerEdgeAlpha', 0); hold on;
    
    
    data_area1 = exp_data{ind1}.area;
    data_exp1 = exp_data{ind1}.exp;    

    if (dataset == 1)    
        scatter(data_area1, data_exp1, 6, colorB1, 'filled', 'MarkerFaceAlpha', alpha_val(type), 'MarkerEdgeAlpha', 0); hold on;    
    else
        scatter(data_area1, data_exp1, 6, colorG1, 'filled', 'MarkerFaceAlpha', alpha_val(type), 'MarkerEdgeAlpha', 0); hold on;
    end
    
    % (2) Plot prediction from model 
    
    if (dataset == 1)
                
        plot(area_axis, mean_traj{ind1}, '-', 'color', colorB2, 'linewidth', 2); hold on;
        
        plot(exp_data_binned{ind1}.xval, exp_data_binned{ind1}.yval, 'sq', 'color', colorB2, 'linewidth',1.5); hold on;

    else
       
        plot(area_axis, mean_traj{ind1}, '-', 'color', colorG2, 'linewidth', 2); hold on;
        
        plot(exp_data_binned{ind1}.xval, exp_data_binned{ind1}.yval, 'sq', 'color', colorG2, 'linewidth',1.5); hold on;        
            
    end
    
    plot(area_axis, mean_traj{ind2}, '-', 'color', colorY2, 'linewidth', 2); hold on;
    
    plot(exp_data_binned{ind2}.xval, exp_data_binned{ind2}.yval, 'sq', 'color', colorY2,'linewidth', 1.5); hold on;  
    
    xlim([1 10]);
    
    
    cd(save_path);
    %saveas( h3B, strcat(save_name, '_h3B.tif') );
    %saveas( h3B, strcat(save_name, '_h3B.fig') );

    
end

%=========================================================================

function [mean_traj, area_axis] = get_mean_trajectory(simudata)

% Averaging trajectories for repeated optimization results.

num_repeat = length(simudata);

rec = {};  % for six different curves

for type = 1:6
    rec{type} = [];
end

for j = 1:num_repeat
        
    temp = {};
    temp{1} = simudata{j}.GR.F.nu;  
    temp{2} = simudata{j}.GR.G.nu;
    temp{3} = simudata{j}.AF.F.bin_AF_RNAP;
    temp{4} = simudata{j}.AF.G.bin_AF_RNAP;    
    temp{5} = simudata{j}.AF.F.bin_AF_ribo;
    temp{6} = simudata{j}.AF.G.bin_AF_ribo;
    
    for type = 1:6        
        rec{type} = [rec{type}; temp{type}(:,2)'];        
    end
    
end

mean_traj = {};

for type = 1:6    
    mean_traj{type} = mean(rec{type},1);    
end

area_axis = simudata{1}.AF.area_axis;

end

%=========================================================================
