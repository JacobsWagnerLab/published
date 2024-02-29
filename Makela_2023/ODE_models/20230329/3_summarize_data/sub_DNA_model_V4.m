
function [output] = sub_DNA_model_V4(param_temp)

% DNA limitation model

% =========================================================================
%%% (1) Model parameters

%param = z;

r1 = param_temp.r1;
r2 = param_temp.r2;
K1 = param_temp.K1;      % genome / um^3
K2 = param_temp.K2;      % mRNA / um^3
d = param_temp.d;
Zs = param_temp.Zs;      % DNA mean concentration (genome / um^3)
c = param_temp.c;       
cp = param_temp.cp;
xini = param_temp.xini;
yini = param_temp.yini;


%%% (2) Models:

% (2-1) DNA scales with cell volume and protein

p = ( r1* Zs ) / (K1 + Zs);
m = K2 * c;

F1 = @(x) p* x(2) - d* x(1);
F2 = @(x) ( r2 * x(2) ) / ( m*(x(2)/x(1)) + 1 ) ;

F = @(x) [F1(x) F2(x)]';
Ft = @(t,x) F(x);

% (2-2) DNA content stay constant

Z0 = 1;
a = (K1 * c) / Z0;
b = K2 * c;

G1 = @(x) ( ( r1 * x(2) ) / ( 1 + a*x(2) ) ) -  d* x(1);
G2 = @(x) ( r2 * x(2) ) / ( b*(x(2)/x(1)) + 1 ) ;

G = @(x) [G1(x) G2(x)]';
Gt = @(t,x) G(x);

% =========================================================================

ini = [xini yini]';
Tstep = 1;
TmaxA = 2500;
TmaxB = 2500;
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

[t_solF, x_solF] = ode23(Ft, [0:Tstep:TmaxA], ini, opt); 
dataF = analyze_traj(t_solF, x_solF, ini, Tstep, param_temp);

dataF.phiP = yini/xini;
dataF.phi = dataF.y(end)/dataF.x(end);           % balanced fraction
dataF.lambdaB = p* dataF.phi - d;                % theoretical formula 
dataF.tau = log(2)/dataF.lambdaB;

% (II) Calculate DNA limited case; start with balanced fraction

ini = yini*[1/dataF.phi 1]';

[t_solG, x_solG] = ode23(Gt, [0:Tstep:TmaxB], ini, opt); 
dataG = analyze_traj(t_solG, x_solG, ini, Tstep, param_temp);

% ----------------------------------------------------------------------

output = {};
output.F = dataF;
output.G = dataG;
output.param = param_temp;
output.paramB = paramB;


end

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
