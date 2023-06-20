
function [data_summaryC] = optimized_simulation(dataset_flag, data_summaryB)

% Using optimized parameter to perform simulation
% Manaully load "data_summaryB.mat"

%dataset_flag = 3;

param_mat = data_summaryB.param_summaryB.mat;

simudata = {};
repeat_num = size(param_mat,2);

for m = 1:repeat_num

    param_temp = param_mat(:,m);
    param_opt = load_optimized_param(param_temp);
    simudata{m} = sub_simulate_single_param(param_opt, dataset_flag);
    
    simudata{m}.GR = get_GR_from_traj(simudata{m});
    simudata{m}.AF = get_AF_from_traj(simudata{m});    
    
end

% Save the results by appending to the data_summaryB and call it
% data_summaryC

data_summaryC = data_summaryB;
data_summaryC.simudata = simudata;

%save_dir = '/Users/wei-hsiang/Desktop/DNA content model/post_analysis/3_post_analysis_C/';
%cd(save_dir); save('data_summaryC', 'data_summaryC');



%

cB = 0.5*[1 1 1];
cY = [0.8 0.7 0];

figure;

for m = 1:repeat_num
    
    subplot(6,5,m);
    
    area = simudata{m}.AF.area_axis;
    ribo_F = simudata{m}.AF.F.bin_AF_ribo;
    ribo_G = simudata{m}.AF.G.bin_AF_ribo;
    
    plot(area, ribo_F(:,2), '-', 'color', cB);  hold on;
    plot(area, ribo_G(:,2), '-', 'color', cY);
    
    ylim([0 1]);
    
    
end



figure;

for m = 1:repeat_num
    
    %subplot(6,5,m);
    
    area = simudata{m}.AF.area_axis;
    ribo_F = simudata{m}.AF.F.bin_AF_ribo;
    ribo_G = simudata{m}.AF.G.bin_AF_ribo;
    
    plot(area, ribo_F(:,2), '-', 'color', cB);  hold on;
    plot(area, ribo_G(:,2), '-', 'color', cY);  hold on;
    
    ylim([0 1]);    
    
end

end


%========================================================================

function [output] = get_GR_from_traj(dataS)

output = {};

F = dataS.traj.F;
G = dataS.traj.G;

area_axis = 0.5*(1:30)';  % um^2

output.F.nu = [area_axis interp1(F.A2, F.Anu, area_axis)];
output.G.nu = [area_axis interp1(G.A2, G.Anu, area_axis)];

end

%========================================================================

function [output] = get_AF_from_traj(dataS)

output = {};

F = dataS.traj.F;
G = dataS.traj.G;
param = dataS.traj.param;

DNA_conc = param.Zs;   % genome copy / um^3
DNA_num = 1;           % genome copy / cell

output.F.area = F.A;

output.F.act_RNAP = ones(length(F.t),1) * ( DNA_conc / (DNA_conc + param.K1 ) );

output.F.act_ribo = F.x ./( param.K2 * param.c * F.y + F.x  );

output.G.area = G.A;

output.G.act_RNAP = DNA_num ./ ( DNA_num + (param.K1 * param.c * G.y) );

output.G.act_ribo = G.x ./( param.K2 * param.c * G.y + G.x );

area_axis = 0.5*(1:30)';  % um^2

output.F.bin_AF_RNAP = [area_axis interp1(output.F.area, output.F.act_RNAP, area_axis)];
output.G.bin_AF_RNAP = [area_axis interp1(output.G.area, output.G.act_RNAP, area_axis)];
output.F.bin_AF_ribo = [area_axis interp1(output.F.area, output.F.act_ribo, area_axis)];
output.G.bin_AF_ribo = [area_axis interp1(output.G.area, output.G.act_ribo, area_axis)];
output.area_axis = area_axis;

end
