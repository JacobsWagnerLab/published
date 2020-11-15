
% SRN-equivalent model for cross-feeding population in chemostat
%
% This code implement the SRN-equivalent model for cross-feeding
% population in chemostat. In equation, x(1), x(2), x(3) are three
% species of microorganisms, x(4), x(5), x(6) are the cross-feeding
% metabolites they produce, and x(7) are the limited nutreint in chemostate
% (denoted by S in the SI, Method section). The quantities of x(t) here
% denotes the total mass of cell or metabolite in a chemostate with 
% unit volume (V=1). 
%
% By simulating this system with different dilution rate D's, 
% we can generate different dynamics as shown in Fig 4E. 
%
% WHLin, 2020 Aug


DR = 0.5;   % Dilution rate 

am = [0.383 0.660 0.235; 0.351 0.354 0.667; 0.665 0.024 0.270]; % competition matrix
cm = [0.084 0.003 0.049; 0.053 0.032 0.063; 0.063 0.022 0.019]; % cross-feeding matrix
xini = [0.018 0.012 0.038 0.592 0.296 0.044 0.5];               % initial condition

param.am = am;
param.cm = cm;

rec = {};
[t_sol, x_sol] = confirm_plot(xini, param, DR);

figure; 
plot(t_sol(8000:10000), x_sol(8000:10000, 1:6), '-');

% ============================================================== %

function [t_sol, x_sol] = confirm_plot(xini, param, DR)

a = param.am;
c = param.cm;

Tstep = 1;
Tmax = 10000;

% External nutrient influx param
r_max = 2;    % maximal uptake rate   (1/time)
KS = 1;       % uptake rate parameter (g/L)  
S_in = 1.5;   % total nutreint influx (g/time) 

% Cross-feeding prodution
b = 4*[0.1 0.1 0.1];    

% Corss-feeding influx param
Km = 0.1;        

% ================================================== %

Jin = @(x) r_max* ( x(7) / (x(7) + KS) );
Jr1 = @(x) Jin(x) * x(1);
Jr2 = @(x) Jin(x) * x(2);
Jr3 = @(x) Jin(x) * x(3);

xv = @(x) [x(1) x(2) x(3)]';

J1 = @(x) Jr1(x) * ( 1 - a(1,:) * xv(x) ) ;
J2 = @(x) Jr2(x) * ( 1 - a(2,:) * xv(x) ) ;
J3 = @(x) Jr3(x) * ( 1 - a(3,:) * xv(x) ) ;

% Secretory fluxes

J4 = @(x) b(1) * x(1);
J5 = @(x) b(2) * x(2);
J6 = @(x) b(3) * x(3);

% Cross-feeding flux

W12 = @(x) c(1,2) * x(1) * ( x(5)/(Km + x(5)) );
W13 = @(x) c(1,3) * x(1) * ( x(6)/(Km + x(6)) );

W21 = @(x) c(2,1) * x(2) * ( x(4)/(Km + x(4)) );
W23 = @(x) c(2,3) * x(2) * ( x(6)/(Km + x(6)) );

W31 = @(x) c(3,1) * x(3) * ( x(4)/(Km + x(4)) );
W32 = @(x) c(3,2) * x(3) * ( x(5)/(Km + x(5)) );

%

F1 = @(x) J1(x) - J4(x) + W12(x) + W13(x) - DR* x(1);
F2 = @(x) J2(x) - J5(x) + W21(x) + W23(x) - DR* x(2);
F3 = @(x) J3(x) - J6(x) + W31(x) + W32(x) - DR* x(3);

F4 = @(x) J4(x) - W21(x) - W31(x) - DR* x(4);
F5 = @(x) J5(x) - W12(x) - W32(x) - DR* x(5);
F6 = @(x) J6(x) - W13(x) - W23(x) - DR* x(6);

F7 = @(x) DR* S_in -  J1(x) - J2(x) - J3(x) - DR* x(7);

Ft = @(t,x) [F1(x) F2(x) F3(x) F4(x) F5(x) F6(x) F7(x)]';

optA = odeset('RelTol', 10^(-4),'AbsTol', 10^(-4));

[t_sol, x_sol] = ode45(Ft,[0:Tstep:Tmax], xini, optA);


end

