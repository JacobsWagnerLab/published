

dataS = simuPQv2;

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

figure;
plot(output.F.area, output.F.act_RNAP, '.');

output.F.bin_AF_RNAP = [area_axis interp1(output.F.area, output.F.act_RNAP, area_axis)];
output.G.bin_AF_RNAP = [area_axis interp1(output.G.area, output.G.act_RNAP, area_axis)];
output.F.bin_AF_ribo = [area_axis interp1(output.F.area, output.F.act_ribo, area_axis)];
output.G.bin_AF_ribo = [area_axis interp1(output.G.area, output.G.act_ribo, area_axis)];
output.area_axis = area_axis;
