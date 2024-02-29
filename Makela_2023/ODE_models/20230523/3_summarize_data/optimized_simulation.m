
function [data_summaryC] = optimized_simulation(dataset_flag, data_summaryB)

% Using optimized parameter to perform simulation

param_mat = data_summaryB.param_summaryB.mat;

simudata = {};
repeat_num = size(param_mat,2);

for m = 1:repeat_num

    param_temp = param_mat(:,m);
    param_opt = load_optimized_param(param_temp);
    simudata{m} = sub_simulate_single_param(param_opt, dataset_flag);
    
    simudata{m}.GR = get_GR_from_traj(simudata{m});
    simudata{m}.AF = get_AF_from_traj_PQv2(simudata{m});    
    
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
% 
% function [output] = get_AF_from_traj(dataS)
% 
% output = {};
% 
% F = dataS.traj.F;
% G = dataS.traj.G;
% param = dataS.traj.param;
% 
% DNA_conc = param.Zs;   % genome copy / um^3
% DNA_num = 1;           % genome copy / cell
% 
% output.F.area = F.A;
% 
% output.F.act_RNAP = ones(length(F.t),1) * ( DNA_conc / (DNA_conc + param.K1 ) );
% 
% output.F.act_ribo = F.x ./( param.K2 * param.c * F.y + F.x  );
% 
% output.G.area = G.A;
% 
% output.G.act_RNAP = DNA_num ./ ( DNA_num + (param.K1 * param.c * G.y) );
% 
% output.G.act_ribo = G.x ./( param.K2 * param.c * G.y + G.x );
% 
% area_axis = 0.5*(1:30)';  % um^2
% 
% output.F.bin_AF_RNAP = [area_axis interp1(output.F.area, output.F.act_RNAP, area_axis)];
% output.G.bin_AF_RNAP = [area_axis interp1(output.G.area, output.G.act_RNAP, area_axis)];
% output.F.bin_AF_ribo = [area_axis interp1(output.F.area, output.F.act_ribo, area_axis)];
% output.G.bin_AF_ribo = [area_axis interp1(output.G.area, output.G.act_ribo, area_axis)];
% output.area_axis = area_axis;
% 
% end

%========================================================================

function [output] = get_AF_from_traj_PQv2_bk(dataS)

output = {};

F = dataS.traj.F;
G = dataS.traj.G;
param = dataS.traj.param;

output.F.area = F.A;

output.F.act_RNAP = 0.5* (F.A>0);
output.F.act_ribo = F.x ./( param.K2 * param.c * F.y + F.x  );

output.G.area = G.A;

% (*) Define active fraction using PQ model. Values for M9GlyCAA
% k1 = 7.5* 10^4;  % 1 / (M.sec)
% km1 = 1.1;       % 1 / sec
% k2 = 0.012;      % 1 / sec
% k3 = 0.0083;     % 1 / sec

nQ = 3000;
mRNAP = 1.05* 10^(-10); 
theta_RNAP = 1.08 * 10^(-3);

Z0 = 1;          % genome copy / cell
Y0 = 6.4* 10^6;
% beta_mRNA = 0.18;
% T_RNAP = 0.33;
% m = K2 * c;

P = (1/param.c)* (theta_RNAP + mRNAP* (G.y - Y0));  % absolute #
Q = ( Z0 * nQ ) ./ ( param.c * G.y );               % absolute #
V = param.c * G.y;      % um^3
Vp = V*10^(-15);       % L
Nav = 6* 10^23;

Pconc = (P/Nav)./Vp;    % mole/V = M
Qconc = (Q/Nav)./Vp;    % mole/V = M

A = (param.k1/(param.km1 + param.k2)) * (1 + (param.k2 / param.k3));
B = (param.k1/(param.km1 + param.k2));

h1 = (1 + A*Qconc - B*Pconc) ./ (2*B*Pconc);
h2 = 1 ./(B*Pconc);

output.G.act_RNAP = 1 - ( sqrt( (h1.^2) + h2 ) - h1 );

output.G.act_ribo = G.x ./( param.K2 * param.c * G.y + G.x );

area_axis = 0.5*(1:30)';  % um^2

output.F.bin_AF_RNAP = [area_axis interp1(output.F.area, output.F.act_RNAP, area_axis)];
output.G.bin_AF_RNAP = [area_axis interp1(output.G.area, output.G.act_RNAP, area_axis)];
output.F.bin_AF_ribo = [area_axis interp1(output.F.area, output.F.act_ribo, area_axis)];
output.G.bin_AF_ribo = [area_axis interp1(output.G.area, output.G.act_ribo, area_axis)];
output.area_axis = area_axis;


end

% ========================================================================

function [output] = get_AF_from_traj_PQv2(dataS)
output = {};

F = dataS.traj.F;
G = dataS.traj.G;
param = dataS.traj.param;

nQ = 3000;
b = 430; 
theta_RNAP0 = 1.08 * 10^(-3);
A = (param.k1/(param.km1 + param.k2)) * (1 + (param.k2 / param.k3));
B = (param.k1/(param.km1 + param.k2));

Zbal = 1.4;     % genome copy / um^3
Z0 = 1;          % genome copy / cell
Y0 = 6.4* 10^6;
Nav = 6* 10^23;

% ---------------------------------------------------------------

output.F.area = F.A;

PF = (1/param.c)* theta_RNAP0* (F.x>0);    % #/um^3
QF = Zbal*nQ *(F.x>0) ;                    % #/um^3

PconcF = (PF/Nav) * 10^15;    % mole/V = M
QconcF = (QF/Nav) * 10^15;    % mole/V = M

h1F = (1 + A*QconcF - B*PconcF) ./ (2*B*PconcF);
h2F = 1 ./(B*PconcF);

output.F.act_RNAP = 1 - ( sqrt( (h1F.^2) + h2F ) - h1F );
output.F.act_ribo = F.x ./( param.K2 * param.c * F.y + F.x  );

% ---------------------------------------------------------------

output.G.area = G.A;

PG = (1/param.c)* (theta_RNAP0 + b*param.c*param.cp* (G.y - Y0));  % #/um^3
QG = ( Z0 * nQ ) ./ ( param.c * G.y ) ;                            % #/um^3

PconcG = (PG/Nav) * 10^15;    % mole/V = M
QconcG = (QG/Nav) * 10^15;    % mole/V = M

h1G = (1 + A*QconcG - B*PconcG) ./ (2*B*PconcG);
h2G = 1 ./(B*PconcG);

output.G.act_RNAP = 1 - ( sqrt( (h1G.^2) + h2G ) - h1G );

output.G.act_ribo = G.x ./( param.K2 * param.c * G.y + G.x  );

% ---------------------------------------------------------------

area_axis = 0.5*(1:30)';  % um^2

%figure;plot(output.F.area, output.F.act_RNAP, '.');

output.F.bin_AF_RNAP = [area_axis interp1(output.F.area, output.F.act_RNAP, area_axis)];
output.G.bin_AF_RNAP = [area_axis interp1(output.G.area, output.G.act_RNAP, area_axis)];
output.F.bin_AF_ribo = [area_axis interp1(output.F.area, output.F.act_ribo, area_axis)];
output.G.bin_AF_ribo = [area_axis interp1(output.G.area, output.G.act_ribo, area_axis)];
output.area_axis = area_axis;

end

