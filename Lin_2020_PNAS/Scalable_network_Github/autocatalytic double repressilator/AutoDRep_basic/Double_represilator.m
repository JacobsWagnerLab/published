
% This code simulates autocatalytic double represillator with seven nodes. 
% For equations and network diagram, see Figure 2F and Supplementary Method.
% For integration scheme, see 'ode45' function for description.
% WHLin, 2020 Jan

% ======== Simulation parameters ========

% Double repressilator
% Simulating a growing double repressilator

Tstep = 0.05;  % step size of time points
Tmax = 1000;   % maximal time step
optA = odeset('RelTol', 10^(-8), 'AbsTol', 10^(-8));  % error tolerance

% Note: the Tstep here is NOT the real integration step size,
% see 'ode45' Matlab document for how integration step size is 
% determined from error tolerance. 

y0 = [0.21 0.2 0.1 0.1 0.2 0.1 0.1];   % Initial condition
y0 = y0/sum(y0);            % Normalized initial condition   

% ======== Model parameters ========

b = [NaN 0.9 0.7 1.0 0.6 0.4 1.2];
c = 10;
d = [NaN 0.2 0.5 0.3 0.1 0.4 0.6];

Ma = 100;
Mb = 100;

Ka = 1500;
alpha = 500;
Kb = alpha*Ka;

theta_a = 4;
theta_b = 4;
phi_a = 4;
phi_b = 4;

% ======== Model description ========

% Repressilator A: node {2,3,4}
% Repressilator B: node {5,6,7}

Iab1 = @(y) 1 / ( 1 + Ma* ( y(2)^phi_a ) );  % Inhibition from repressilator A to repressilator B
Iba1 = @(y) 1 / ( 1 + Mb* ( y(5)^phi_b ) );  % Inhibition from repressilator B to repressilator A

Jin = @(y) b(2)*y(2) + b(3)*y(3) + b(4)*y(4) + b(5)*y(5) + b(6)*y(6) + b(7)*y(7);
J1 = @(y) Jin(y);

J2 = @(y) Iba1(y)* c* y(1) / (1 + Ka* ( y(3)^theta_a ) ) ;
J3 = @(y)          c* y(1) / (1 + Ka* ( y(4)^theta_a ) ) ;
J4 = @(y)          c* y(1) / (1 + Ka* ( y(2)^theta_a ) ) ;

J5 = @(y) Iab1(y)* c* y(1) / (1 + Kb* ( y(6)^theta_b ) ) ;
J6 = @(y)          c* y(1) / (1 + Kb* ( y(7)^theta_b ) ) ;
J7 = @(y)          c* y(1) / (1 + Kb* ( y(5)^theta_b ) ) ;

Jd2 = @(y) d(2)*y(2);
Jd3 = @(y) d(3)*y(3);
Jd4 = @(y) d(4)*y(4);
Jd5 = @(y) d(5)*y(5);
Jd6 = @(y) d(6)*y(6);
Jd7 = @(y) d(7)*y(7);

F1 = @(y) J1(y) - J2(y) - J3(y) - J4(y) - J5(y) - J6(y) - J7(y);
F2 = @(y) J2(y) - Jd2(y);
F3 = @(y) J3(y) - Jd3(y);
F4 = @(y) J4(y) - Jd4(y);
F5 = @(y) J5(y) - Jd5(y);
F6 = @(y) J6(y) - Jd6(y);
F7 = @(y) J7(y) - Jd7(y);

F = @(y) [F1(y) F2(y) F3(y) F4(y) F5(y) F6(y) F7(y)]';
mu = @(y) sum(F(y));

% Project to the simplex space
G = @(t,y) F(y) - mu(y).*y; 

% Integrating the model
[t_sol, y_sol] = ode45(G,[0:Tstep:Tmax], y0, optA);


% =========  Visualizing solution trajectory =======

time_rg = ceil( [0.7*Tmax Tmax]/Tstep ) ;  % fraction of trajecotry to be plot

t_eq = t_sol(time_rg(1)+1:time_rg(2));     % time points for plotting
y_eq = y_sol(time_rg(1)+1:time_rg(2),:);   % trajectory Y(t) for plitting


h1 = figure('unit', 'normalized', 'outerposition', [0 0 0.8 0.8]);

subplot(2,4,[1 4]);  
plot(t_eq, y_eq(:,1), t_eq, y_eq(:,2), t_eq, y_eq(:,5));
xlabel('time');  ylabel('y(t)');

subplot(2,4,[5 6]); 
plot(y_eq(:,2), y_eq(:,3), '-');
xlim([0 0.4]);  
ylim([0 0.4]);
xlabel('y_2');  ylabel('y_3');

subplot(2,4,[7 8]); 
plot(y_eq(:,2), y_eq(:,5), '-');
xlim([0 0.4]);  
ylim([0 0.4]);
xlabel('y_2');  ylabel('y_5');
