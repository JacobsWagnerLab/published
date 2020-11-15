
% This script siumlation autocatalytic repressilator with four nodes.
% The simulation generates plot Figure 2B, 2C in the main text. 
% See Supplementary Method for the equations in the model.
% WHLin, 2020 Jan.

% ============ Simulation Parameters ============ %

Tstep = 0.1;  % time step in the trajectory
Tmax = 40;    % maximal simulation time
optA = odeset('RelTol', 10^(-6),'AbsTol', 10^(-6));   % error tolerance

% =============== Model Parameters ============== %

b = [NaN 0.5 0.4 0.3];
csyn = 20;
d = [NaN 0.25 0.3 0.2];

theta = 1;   % Hill Coefficient
K = 500;     % Repression strength

% =============== Model Description ============== %

% (1) Original space X(t)

Jin = @(x) b(2)*x(2) + b(3)*x(3) + b(4)*x(4);
J1 = @(x) Jin(x);
N = @(x) sum(x);

J2 = @(x) csyn * x(1) / (1 + K* ( ( x(3)/ N(x) ) ^theta ) ) ;
J3 = @(x) csyn * x(1) / (1 + K* ( ( x(4)/ N(x) ) ^theta ) ) ;
J4 = @(x) csyn * x(1) / (1 + K* ( ( x(2)/ N(x) ) ^theta ) ) ;

Jd2 = @(x) d(2)*x(2);
Jd3 = @(x) d(3)*x(3);
Jd4 = @(x) d(4)*x(4);

F1x = @(x) J1(x) - J2(x) - J3(x) - J4(x);
F2x = @(x) J2(x) - Jd2(x);
F3x = @(x) J3(x) - Jd3(x);
F4x = @(x) J4(x) - Jd4(x);

Fori = @(x) [F1x(x) F2x(x) F3x(x) F4x(x)]';

Ft = @(t,x) Fori(x);

% (2) Rescaled space Y(t)

Jin = @(y) b(2)*y(2) + b(3)*y(3) + b(4)*y(4);
J1 = @(y) Jin(y);

J2 = @(y) csyn * y(1) / (1 + K* ( y(3)^theta ) ) ;
J3 = @(y) csyn * y(1) / (1 + K* ( y(4)^theta ) ) ;
J4 = @(y) csyn * y(1) / (1 + K* ( y(2)^theta ) ) ;

Jd2 = @(y) d(2)*y(2);
Jd3 = @(y) d(3)*y(3);
Jd4 = @(y) d(4)*y(4);

F1 = @(y) J1(y) - J2(y) - J3(y) - J4(y);
F2 = @(y) J2(y) - Jd2(y);
F3 = @(y) J3(y) - Jd3(y);
F4 = @(y) J4(y) - Jd4(y);

F = @(y) [F1(y) F2(y) F3(y) F4(y)]';
mu = @(y) sum(F(y));

Gt = @(t,y) F(y) - mu(y).*y; % Project to the simplex space

 
% Simulating the equation 
x0 = [0.1 0.3 0.3 0.3];  % initial condition
y0 = x0/sum(x0);         % normalized initial condition

[t_solx, x_sol] = ode23(Ft,[0:Tstep:Tmax], x0, optA);
[t_soly, y_sol] = ode45(Gt,[0:Tstep:Tmax], y0, optA);

% =============== Model Description ============== %

h = figure('position', [1 1 300 500]);

subplot(211); 
semilogy(t_solx, x_sol, '-', 'LineWidth', 1.5);
ylabel('biomass vector x(t)');
xlabel('time'); 

subplot(212); 
plot(t_soly, y_sol, '-', 'LineWidth', 1.5);
ylabel('normalized biomass vector y(t)');
xlabel('time'); 

% End of the script

