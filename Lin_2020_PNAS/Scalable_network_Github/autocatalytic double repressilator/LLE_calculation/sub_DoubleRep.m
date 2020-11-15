
% This code simulates autocatalytic double represillator with seven nodes. 
% For equations and network diagram, see Figure 2F and Supplementary Method.
% For integration scheme, see 'ode45' function for description.
% This model include scalable noise, allow setting noise level.
% This is a modified subfunction for simulating various parameter (alpha)

% WHLin, 2020 Jan

% ======== Input parameters ========
% - y0: initial condition
% - Tstep: intervals of time points for the simulated trajectory 
% Note that this is not the step size used for integration
% see Matlab ode45() for description

% - Tmax: maximal simulation time
% - param: specified parameter, (alpha)
% - noise_level: level of scalable noise in the system

% Example for inpur parameters:
% y0 = [10 0.2 0.15 0.15 0.2 0.15 0.15];
% y0 = y0/sum(y0);
% Tstep = 0.1;
% Tmax = 500;
% param = 100;
% noise_level = 0.05;

% ======== Output results ========
% - t_sol: time points solution trajecotry
% - y_sol: solution trajectory Y(t)
% - GR_rec: long-term growth rate calculated from Y(t)

% ===================================================================== %

function [t_sol, y_sol, GR_rec] = sub_DoubleRep(y0, Tstep, Tmax, param, noise_level)

% ======== Model parameters ========

b = 2*[0 0.45 0.35 0.5 0.3 0.2 0.6];
c = 10*[0 1 1 1 1 1 1];
d = 1*[0 0.2 0.5 0.3 0.1 0.4 0.6];

Ma = 100;
Mb = 100;

Ka = 1500;
Kb = param*Ka;

theta_a = 4;
theta_b = 4;
phi_a = 4;
phi_b = 4;

% ======== Model description ========

Iab1 = @(y) 1 / ( 1 + Ma* ( y(2)^phi_a ) );
Iba1 = @(y) 1 / ( 1 + Mb* ( y(5)^phi_b ) );

Jin = @(y) b(2)*y(2) + b(3)*y(3) + b(4)*y(4) + b(5)*y(5) + b(6)*y(6) + b(7)*y(7);
J1 = @(y) Jin(y);

J2 = @(y) Iba1(y)*c(2)* y(1) / (1 + Ka* ( y(3)^theta_a ) ) ;
J3 = @(y)         c(3)* y(1) / (1 + Ka* ( y(4)^theta_a ) ) ;
J4 = @(y)         c(4)* y(1) / (1 + Ka* ( y(2)^theta_a ) ) ;

J5 = @(y) Iab1(y)*c(5)* y(1) / (1 + Kb* ( y(6)^theta_b ) ) ;
J6 = @(y)         c(6)* y(1) / (1 + Kb* ( y(7)^theta_b ) ) ;
J7 = @(y)         c(7)* y(1) / (1 + Kb* ( y(5)^theta_b ) ) ;

Jd2 = @(y) d(2)*y(2);
Jd3 = @(y) d(3)*y(3);
Jd4 = @(y) d(4)*y(4);
Jd5 = @(y) d(5)*y(5);
Jd6 = @(y) d(6)*y(6);
Jd7 = @(y) d(7)*y(7);

Z1 = @(y) noise_level* y(1)* rand(1,1);
Z2 = @(y) noise_level* y(2)* rand(1,1);
Z3 = @(y) noise_level* y(3)* rand(1,1);
Z4 = @(y) noise_level* y(4)* rand(1,1);
Z5 = @(y) noise_level* y(5)* rand(1,1);
Z6 = @(y) noise_level* y(6)* rand(1,1);
Z7 = @(y) noise_level* y(7)* rand(1,1);

F1 = @(y) Z1(y) + J1(y) - J2(y) - J3(y) - J4(y) - J5(y) - J6(y) - J7(y);
F2 = @(y) Z2(y) + J2(y) - Jd2(y) ;
F3 = @(y) Z3(y) + J3(y) - Jd3(y);
F4 = @(y) Z4(y) + J4(y) - Jd4(y);
F5 = @(y) Z5(y) + J5(y) - Jd5(y);
F6 = @(y) Z6(y) + J6(y) - Jd6(y);
F7 = @(y) Z7(y) + J7(y) - Jd7(y);

F = @(y) [F1(y) F2(y) F3(y) F4(y) F5(y) F6(y) F7(y)]';

mu = @(y) J1(y) - Jd2(y) - Jd3(y) - Jd4(y) - Jd5(y) - Jd6(y) - Jd7(y);
nu = @(y) sum(F(y));

% Project to the simplex space
G = @(t,y) F(y) - nu(y).*y;

[t_sol, y_sol] = ode45(G,[0:Tstep:Tmax], y0);


% ======== Calculating long-term growth rate ========

L = size(y_sol,1);
GR_rec = NaN(L,1);

for j = 1:L
    
    GR_rec(j) = mu(y_sol(j,:));

end

% End of the script


