
function [output] = sub_DNA_model_PQv2(param_temp)
%param_temp = paramI1;

% DNA limitation model

% =========================================================================
%%% (1) Model parameters

r1 = param_temp.r1;
r2 = param_temp.r2;
K1 = param_temp.K1;      % genome / um^3
K2 = param_temp.K2;      % genome / um^2
d = param_temp.d;
Zs = param_temp.Zs;      % DNA mean concentration (genome / um^3)
c = param_temp.c;       
cp = param_temp.cp;
xini = param_temp.xini;
yini = param_temp.yini;

k1 = param_temp.k1;
km1 = param_temp.km1;
k2 = param_temp.k2;
k3 = param_temp.k3;

%%% (2) Models:

% (2-0) Define active fraction using PQ model. Values for M9GlyCAA

nQ = 3000;
b = 430; 
theta_RNAP0 = 1.08 * 10^(-3);
A = (k1/(km1+k2)) * (1 + (k2/k3));
B = (k1/(km1+k2));

Zbal = 1.4;   % genome#/um^3
Z0 = 1;       % genome#/cell

Y0 = 6.4* 10^6;
beta_mRNA = 0.18;
T_RNAP = 0.33;
Nav = 6* 10^23;
m = K2 * c;

%-----------------------------------------------------------------
% (1) Balanced growth

r1F = @(x) theta_RNAP0 * (beta_mRNA/T_RNAP);

PF = @(x) (1/c)* theta_RNAP0;   % #/um^3
QF = @(x) Zbal*nQ ;             % #/um^3

PconcF = @(x) (PF(x)/Nav) * 10^15;    % mole/L = M
QconcF = @(x) (QF(x)/Nav) * 10^15;    % mole/L = M

h1F = @(x) (1 + A*QconcF(x) - B*PconcF(x)) / (2*B*PconcF(x));
h2F = @(x) 1 /(B*PconcF(x));

AF_RNAP_F = @(x) 1 - ( sqrt( h1F(x)*h1F(x) + h2F(x) ) - h1F(x) );

F1 = @(x) r1F(x)* AF_RNAP_F(x) * x(2) - d* x(1);
F2 = @(x) ( r2 * x(2) ) / ( m*(x(2)/x(1)) + 1 ) ;

F = @(x) [F1(x) F2(x)]';
Ft = @(t,x) F(x);

%-----------------------------------------------------------------
% (2-2) 1N cell

r1G = @(x) ( theta_RNAP0 + b*c*cp*(x(2)-Y0) )* (beta_mRNA/T_RNAP);

PG = @(x) (1/c)* ( theta_RNAP0 + b*c*cp*(x(2)-Y0) );  % #/um^3
QG = @(x) ( Z0 * nQ ) / ( c * x(2) );                 % #/um^3

PconcG = @(x) (PG(x)/Nav) * 10^15;   % mole/V = M
QconcG = @(x) (QG(x)/Nav) * 10^15;   % mole/V = M

h1G = @(x) (1 + A*QconcG(x) - B*PconcG(x)) / (2*B*PconcG(x));
h2G = @(x) 1 /(B*PconcG(x));

AF_RNAP_G = @(x) 1 - ( sqrt( h1G(x)*h1G(x) + h2G(x) ) - h1G(x) );

G1 = @(x) r1G(x)* AF_RNAP_G(x) * x(2) - d* x(1);
G2 = @(x) ( r2 * x(2) ) / ( m*(x(2)/x(1)) + 1 ) ;

G = @(x) [G1(x) G2(x)]';
Gt = @(t,x) G(x);

% =========================================================================

ini = [xini yini]';
Tstep = 1;
TmaxA = 500;
TmaxB = 500;
opt = odeset('RelTol', 10^(-4),'AbsTol', 10^(-5));

paramB = {};
paramB.Tstep = Tstep;
paramB.TmaxA = TmaxA;
paramB.TmaxB = TmaxB;
paramB.opt = opt;


% (Ia) Calculate balanced growth case (pre-set)

[t_solF, x_solF] = ode23(Ft, [0:Tstep:TmaxA], ini, opt); 
dataF = analyze_traj(t_solF, x_solF, ini, Tstep, param_temp);
dataF.phi = dataF.y(end)/dataF.x(end);  % balanced fraction

% (Ib) Calculate balanced growth case (pre-set)

ini = yini*[1/dataF.phi 1]';

[t_solF, x_solF] = ode45(Ft, [0:Tstep:TmaxA], ini, opt); 
dataF = analyze_traj(t_solF, x_solF, ini, Tstep, param_temp);

dataF.phiP = yini/xini;
dataF.phi = dataF.y(end)/dataF.x(end);             % balanced fraction
%dataF.lambdaB = r1* AF_RNAP_F * dataF.phi - d;    % theoretical formula 
%dataF.tau = log(2)/dataF.lambdaB;

% (II) Calculate DNA limited case; start with balanced fraction

ini = yini*[1/dataF.phi 1]';

[t_solG, x_solG] = ode45(Gt, [0:Tstep:TmaxB], ini, opt); 
dataG = analyze_traj(t_solG, x_solG, ini, Tstep, param_temp);

% ----------------------------------------------------------------------

output = {};
output.F = dataF;
output.G = dataG;
output.param = param_temp;
output.paramB = paramB;

end


% figure;
% subplot(121);
% plot(output.G.t, output.G.y);
% subplot(122);
% plot(output.G.t, output.G.x);

% ====================================================================== %

function [data] = analyze_traj(t_sol, x_sol, xini, Tstep, param_temp)

data = {};
data.x = x_sol(:,1);
data.y = x_sol(:,2);
data.t = t_sol;

data.w = data.x ./data.y;
data.N = data.x + data.y;

% Calculate non-normalized instant. growth rate

temp = [];  % [1] dx/dt; [2] dy/dt [3] x [4] y

for j = 1:size(x_sol,1)-1
        
    temp(j,1) = ( x_sol(j+1,1) + x_sol(j,1) ) / 2;
    temp(j,2) = ( x_sol(j+1,2) + x_sol(j,2) ) / 2;
    temp(j,3) = ( x_sol(j+1,1) - x_sol(j,1) ) / Tstep;
    temp(j,4) = ( x_sol(j+1,2) - x_sol(j,2) ) / Tstep;

end

data.x2 = temp(:,1);
data.y2 = temp(:,2);
data.xnu = temp(:,3);
data.ynu = temp(:,4);

% Cell volume and concentration

c = param_temp.c;
cp = param_temp.cp;

data.V = c*data.y;
data.A = cp*data.y;
 
data.xc = data.x ./ data.V;   % concentration of x
data.yc = data.y ./ data.V;   % concentration of y 

data.A2 = cp* data.y2;
data.Anu = cp* data.ynu;

end

% ====================================================================== %
