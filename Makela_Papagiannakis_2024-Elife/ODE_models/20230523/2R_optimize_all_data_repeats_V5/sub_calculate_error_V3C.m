
function [output] = sub_calculate_error_V3C(GR_data, AF_SMT_data, simudata, dataset_flag)

% test fitting with all data points
% Mar 10th. Use maximal error amount six datasets as the error metric

[AF_simu] = get_AF_from_traj_PQv2(simudata);

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

tempS{1} = [simudata.traj.F.A2 simudata.traj.F.Anu];
tempS{2} = [simudata.traj.G.A2 simudata.traj.G.Anu];
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

tempS{1} = [simudata.traj.F.A2 simudata.traj.F.Anu];
tempS{2} = [simudata.traj.G.A2 simudata.traj.G.Anu];
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
    %comparison_rec{m}.abs_err = abs_err;
    %comparison_rec{m}.rel_err = rel_err;
    
    objective_function(m,1) = nanmean( abs(rel_err) );
    objective_function(m,2) = length(rel_err);
    
end

err = {};
err.table = objective_function;
%err.sum = sum(err_weight.*err.table(:,1));
err.max = sum(err.table(:,1));


output = {};
output.rec = comparison_rec;
output.err = err;


%{
figure;
for m = 1:6
    subplot(3,2,m);
    temp = comparison_rec{m};    
    plot( temp.area, temp.exp, '.' ); hold on;
    plot( temp.area, temp.simu, '.' ); hold off; 
end
%}

end


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

function [output] = get_AF_from_traj_PQv2_bk(dataS)

output = {};

F = dataS.traj.F;
G = dataS.traj.G;
param = dataS.traj.param;

nQ = 3000;
mRNAP = 1.05* 10^(-10); 
theta_RNAP = 1.08 * 10^(-3);
A = (param.k1/(param.km1 + param.k2)) * (1 + (param.k2 / param.k3));
B = (param.k1/(param.km1 + param.k2));

Zmean = 1.4;     % genome copy / um3
Z0 = 1;          % genome copy / cell
Y0 = 6.4* 10^6;
Nav = 6* 10^23;

% ---------------------------------------------------------------

output.F.area = F.A;

VF = param.c * F.y;      % um^3
VpF = VF*10^(-15);       % L

Zbal = Zmean * VF;
PF = (1/param.c)* (theta_RNAP + mRNAP* (F.y - Y0));  % absolute #
QF = ( Zbal * nQ ) ./ ( param.c * F.y );               % absolute #

PconcF = (PF/Nav)./VpF;    % mole/V = M
QconcF = (QF/Nav)./VpF;    % mole/V = M

h1F = (1 + A*QconcF - B*PconcF) ./ (2*B*PconcF);
h2F = 1 ./(B*PconcF);

output.F.act_RNAP = 1 - ( sqrt( (h1F.^2) + h2F ) - h1F );
output.F.act_ribo = F.x ./( param.K2 * param.c * F.y + F.x  );

% ---------------------------------------------------------------

output.G.area = G.A;

VG = param.c * G.y;      % um^3
VpG = VG*10^(-15);       % L

PG = (1/param.c)* (theta_RNAP + mRNAP* (G.y - Y0));  % absolute #
QG = ( Z0 * nQ ) ./ ( param.c * G.y );               % absolute #

PconcG = (PG/Nav)./ VpG;    % mole/V = M
QconcG = (QG/Nav)./ VpG;    % mole/V = M

h1G = (1 + A*QconcG - B*PconcG) ./ (2*B*PconcG);
h2G = 1 ./(B*PconcG);

output.G.act_RNAP = 1 - ( sqrt( (h1G.^2) + h2G ) - h1G );

output.G.act_ribo = G.x ./( param.K2 * param.c * G.y + G.x  );

% ---------------------------------------------------------------

area_axis = 0.5*(1:30)';  % um^2

z = [length(output.F.area) length(output.F.act_RNAP)];

figure;
plot(output.F.area, output.F.act_RNAP, '.');
close all;

output.F.bin_AF_RNAP = [area_axis interp1(output.F.area, output.F.act_RNAP, area_axis)];
output.G.bin_AF_RNAP = [area_axis interp1(output.G.area, output.G.act_RNAP, area_axis)];
output.F.bin_AF_ribo = [area_axis interp1(output.F.area, output.F.act_ribo, area_axis)];
output.G.bin_AF_ribo = [area_axis interp1(output.G.area, output.G.act_ribo, area_axis)];
output.area_axis = area_axis;

end

% =======================================================================

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

% =======================================================================
