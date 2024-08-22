
% Compare DNA models

dataset_flag = 1;
paramI = load_literature_param_V5();
paramI1 = paramI{dataset_flag};

simuV4 = sub_DNA_model_V4(paramI1);
simuPQv2 = sub_DNA_model_PQv2(paramI1);

outputV4 = get_AF_from_traj_V4(simuV4);
outputPQ = get_AF_from_traj_PQv2(simuPQv2);


figure;

subplot(131);

tempF = simuV4.F;
tempG = simuV4.G;

plot(tempF.A2, tempF.Anu, 'b-'); hold on;
plot(tempG.A2, tempG.Anu, 'r-'); hold on;

tempFN = simuPQv2.F;
tempGN = simuPQv2.G;

plot(tempFN.A2, tempFN.Anu, 'b-.'); hold on;
plot(tempGN.A2, tempGN.Anu, 'r-.');

xlim([0 15]);


%

subplot(132);

plot(outputV4.area_axis, outputV4.F.bin_AF_RNAP(:,2), 'b-'); hold on;
plot(outputV4.area_axis, outputV4.G.bin_AF_RNAP(:,2), 'r-'); hold on;

plot(outputPQ.area_axis, outputPQ.F.bin_AF_RNAP(:,2), 'bo'); hold on;
plot(outputPQ.area_axis, outputPQ.G.bin_AF_RNAP(:,2), 'ro'); 


subplot(133);

plot(outputV4.area_axis, outputV4.F.bin_AF_ribo(:,2), 'b-'); hold on;
plot(outputV4.area_axis, outputV4.G.bin_AF_ribo(:,2), 'r-'); hold on;

plot(outputPQ.area_axis, outputPQ.F.bin_AF_ribo(:,2), 'bo'); hold on;
plot(outputPQ.area_axis, outputPQ.G.bin_AF_ribo(:,2), 'ro'); 

% 
% subplot(122);
% 
% tempF = simuV4n.F;
% tempG = simuV4n.G;
% 
% plot(tempF.A2, tempF.Anu, 'b-'); hold on;
% plot(tempG.A2, tempG.Anu, 'r-'); hold on;
% 
% tempFN = simuPQv2n.F;
% tempGN = simuPQv2n.G;
% 
% plot(tempFN.A2, tempFN.Anu, 'b-.'); hold on;
% plot(tempGN.A2, tempGN.Anu, 'r-.');
% 
% xlim([0 15]);

%=====================================================================

function [output] = get_AF_from_traj_PQv2_bk(dataS)

output = {};

F = dataS.F;
G = dataS.G;
param = dataS.param;

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
Vp = V* 10^(-15);       % L
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

%=====================================================================

function [output] = get_AF_from_traj_V4(dataS)

output = {};

F = dataS.F;
G = dataS.G;
param = dataS.param;

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

%=====================================================================

function [output] = get_AF_from_traj_PQv2(dataS)

output = {};
F = dataS.F;
G = dataS.G;
param = dataS.param;

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
