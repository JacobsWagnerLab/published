

% test fitting with all data points

dataset_flag = 3;

[param0] = load_literature_param_V2();
[simudata0] = sub_simulate_single_param(param0{dataset_flag}, dataset_flag);
[AF_simu] = get_AF_from_traj(simudata0);

% Compare experimental and simulated data

tempE = {};  % experimental data
tempS = {};  % simulated data

if (dataset_flag == 1)  % M9GlyCAA

ind = dataset_flag*2;    
    
tempE{1} = GR_data{ind-1}.rawdata(:,1:2);  % cell area (A); dA/dt
tempE{2} = GR_data{ind}.rawdata(:,1:2);
tempE{3} = AF_SMT_data{dataset_flag}.ftsZ.RNAP.tableC;  % cell area (A), AF 
tempE{4} = AF_SMT_data{dataset_flag}.oriC.RNAP.tableC;
tempE{5} = AF_SMT_data{dataset_flag}.ftsZ.ribo.tableC;
tempE{6} = AF_SMT_data{dataset_flag}.oriC.ribo.tableC;

tempS{1} = [simudata0.traj.F.A2 simudata0.traj.F.Anu];
tempS{2} = [simudata0.traj.G.A2 simudata0.traj.G.Anu];
tempS{3} = [AF_simu.F.area AF_simu.F.act_RNAP];
tempS{4} = [AF_simu.G.area AF_simu.G.act_RNAP];
tempS{5} = [AF_simu.F.area AF_simu.F.act_ribo];
tempS{6} = [AF_simu.G.area AF_simu.G.act_ribo];

elseif (dataset_flag > 1)  % M9GlyCAA or M9Ala

ind = dataset_flag*2;    
    
tempE{1} = GR_data{ind-1}.rawdata(:,1:2);  % cell area (A); dA/dt
tempE{2} = GR_data{ind}.rawdata(:,1:2);
tempE{3} = AF_SMT_data{dataset_flag}.WT.RNAP.tableC;  % cell area (A), AF 
tempE{4} = AF_SMT_data{dataset_flag}.oriC.RNAP.tableC;
tempE{5} = AF_SMT_data{dataset_flag}.WT.ribo.tableC;
tempE{6} = AF_SMT_data{dataset_flag}.oriC.ribo.tableC;

tempS{1} = [simudata0.traj.F.A2 simudata0.traj.F.Anu];
tempS{2} = [simudata0.traj.G.A2 simudata0.traj.G.Anu];
tempS{3} = [AF_simu.F.area AF_simu.F.act_RNAP];
tempS{4} = [AF_simu.G.area AF_simu.G.act_RNAP];
tempS{5} = [AF_simu.F.area AF_simu.F.act_ribo];
tempS{6} = [AF_simu.G.area AF_simu.G.act_ribo];
    
end

%

comparison_rec = {};
objective_function = NaN(6,2);  % error; statistical N

for m = 1:6

    simu_fit_point = interp1(tempS{m}(:,1), tempS{m}(:,2), tempE{m}(:,1), 'linear', 'extrap');

    abs_err = simu_fit_point - tempE{m}(:,2);
    rel_err = abs_err ./ simu_fit_point;
    
    comparison_rec{m}.area = tempE{m}(:,1);
    comparison_rec{m}.exp = tempE{m}(:,2);
    comparison_rec{m}.simu = simu_fit_point;
    comparison_rec{m}.abs_err = abs_err;
    comparison_rec{m}.rel_err = rel_err;
    
    objective_function(m,1) = nanmean( abs(rel_err) );
    objective_function(m,2) = length(rel_err);
    
end

%{
figure;

for m = 1:6

    subplot(3,2,m);

    temp = comparison_rec{m};    
    plot( temp.area, temp.exp, '.' ); hold on;
    plot( temp.area, temp.simu, '.' ); hold off;
    
end

%}


% =======================================================================

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

% =======================================================================

function [output] = get_AF_from_traj(dataS)

output = {};

F = dataS.traj.F;
G = dataS.traj.G;
param = dataS.traj.param;

DNA_F = param.Zs;
DNA_G = 1;

output.F.area = F.A;

output.F.act_RNAP = ones(length(F.t),1) * ( DNA_F /(DNA_F + param.K1 * param.Vmean) );

output.F.act_ribo = F.x ./(param.K2 * param.c * F.y + F.x);

output.G.area = G.A;

output.G.act_RNAP = DNA_G ./ (param.K1 * param.c * G.y + DNA_G);

output.G.act_ribo = G.x ./(param.K2 * param.c * G.y + G.x);

area_axis = 0.5*(1:30)';  % um^2
%output.F.bin_AF_RNAP = interp1(output.F.area, output.F.act_RNAP, area_axis, 'linear', 'extrap');
%output.G.bin_AF_RNAP = interp1(output.G.area, output.G.act_RNAP, area_axis, 'linear', 'extrap');
%output.F.bin_AF_ribo = interp1(output.F.area, output.F.act_ribo, area_axis, 'linear', 'extrap');
%output.G.bin_AF_ribo = interp1(output.G.area, output.G.act_ribo, area_axis, 'linear', 'extrap');

output.F.bin_AF_RNAP = [area_axis interp1(output.F.area, output.F.act_RNAP, area_axis)];
output.G.bin_AF_RNAP = [area_axis interp1(output.G.area, output.G.act_RNAP, area_axis)];
output.F.bin_AF_ribo = [area_axis interp1(output.F.area, output.F.act_ribo, area_axis)];
output.G.bin_AF_ribo = [area_axis interp1(output.G.area, output.G.act_ribo, area_axis)];
output.area_axis = area_axis;

end

% =======================================================================
