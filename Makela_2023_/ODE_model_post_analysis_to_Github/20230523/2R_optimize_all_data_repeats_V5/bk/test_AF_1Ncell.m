
% Test calculating AF_RNAP

dataS = simudata0;

output = {};

F = dataS.traj.F;
G = dataS.traj.G;
param = dataS.traj.param;

output.F.area = F.A;

output.F.act_RNAP = 0.5* (F.A>0);
output.F.act_ribo = F.x ./( param.K2 * param.c * F.y + F.x  );

output.G.area = G.A;

% (*) Define active fraction using PQ model. Values for M9GlyCAA
k1 = 7.5* 10^4;  % 1 / (M.sec)
km1 = 1.1;       % 1 / sec
k2 = 0.012;      % 1 / sec
k3 = 0.0083;     % 1 / sec
nQ = 3000;
mRNAP = 1.05* 10^(-10); 
theta_RNAP = 1.08 * 10^(-3);

Z0 = 1;          % genome copy / cell
Y0 = 6.4* 10^6;
% beta_mRNA = 0.18;
% T_RNAP = 0.33;
% m = K2 * c;

P = (1/param.c)* (theta_RNAP + mRNAP* (G.y - Y0));
Q = ( Z0 * nQ ) ./ ( param.c * G.y );

A = (k1/(km1+k2)) * (1 + (k2/k3));
B = (k1/(km1+k2));

h1 = (1 + A*Q - B*P) ./ (2*B*P);
h2 = 1 ./(B*P);

output.G.act_RNAP = 1 - ( sqrt( (h1.^2) - h2 ) - h1 );

output.G.act_ribo = G.x ./( param.K2 * param.c * G.y + G.x );

area_axis = 0.5*(1:30)';  % um^2

output.F.bin_AF_RNAP = [area_axis interp1(output.F.area, output.F.act_RNAP, area_axis)];
output.G.bin_AF_RNAP = [area_axis interp1(output.G.area, output.G.act_RNAP, area_axis)];
output.F.bin_AF_ribo = [area_axis interp1(output.F.area, output.F.act_ribo, area_axis)];
output.G.bin_AF_ribo = [area_axis interp1(output.G.area, output.G.act_ribo, area_axis)];
output.area_axis = area_axis;


