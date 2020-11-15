
% This code simulates modified May-Leonard model with three nodes. 
% For equation and network diagram, see FigureS2A and Supplementary Method.
% For integration scheme, see 'ode45' function for description.
% WHLin, 2020 Jan

% The ODE is log-transformed to improve the integration accuracy 
% when component has value closed to 0. This script also includes 
% scalable noise in the model. The result is used in Fig S2E.


% ======== Simulation parameters ========

Tstep = 1;       % step size of time points 
Tmax = 5000;    % maximal time step
optA = odeset('RelTol', 10^(-5),'AbsTol', 10^(-5));  % error tolerance

% Note: the Tstep here is NOT the real integration step size,
% see 'ode45' Matlab document for how integration step size is 
% determined from error tolerance. 

y0 = [2 1 0.5]';    % initial condition
y0 = y0/sum(y0);    % normalized initial condition 

% ======== Model parameters ========

r = [1.5 1.2 0.9];      % maximal growth rate
a = [1.1 1.15 1.12];    % influx coefficients
b = [1 0.9 1.2];        % efflux coefficients

c = [b(1) b(2) b(3); b(3) b(1) b(2); b(2) b(3) b(1)];   % interaction matrix

noise_level = 0.1;      % for deterministic , set noise_level = 0


% ======== Model description ========

% Equation in F(x).  Projected equation in G(y)

N = @(x) x(1)+x(2)+x(3);
y1 = @(x) x(1)/N(x);
y2 = @(x) x(2)/N(x);
y3 = @(x) x(3)/N(x);

Z1 = @(x) noise_level * x(1) * rand(1,1);
Z2 = @(x) noise_level * x(2) * rand(1,1);
Z3 = @(x) noise_level * x(3) * rand(1,1);

F1 = @(x) Z1(x) + ( x(1) * r(1) * ( a(1) - c(1,1)*y1(x) - c(1,2)*y2(x) - c(1,3)*y3(x) ) );
F2 = @(x) Z2(x) + ( x(2) * r(2) * ( a(2) - c(2,1)*y1(x) - c(2,2)*y2(x) - c(2,3)*y3(x) ) );
F3 = @(x) Z3(x) + ( x(3) * r(3) * ( a(3) - c(3,1)*y1(x) - c(3,2)*y2(x) - c(3,3)*y3(x) ) );

mu = @(x) F1(x) + F2(x) + F3(x);
F = @(x) [F1(x) F2(x) F3(x)]';
G = @(y) F(y) - mu(y)*y;
Gt = @(t,x) G(x);

% Log-scale equation for y ( u=log(y) )
% For improving the accuracy when y close to zero, change variable 
% to log-scale for integration

w = @(u) exp(u);
GLog = @(u) (1./w(u)) .* G(w(u));
GLogt = @(t,u) GLog(u);

u0 = log(y0);  % log-transform the initial condition
[t_sol, u_sol] = ode45(GLogt,[0:Tstep:Tmax], u0, optA); 


% ======== Data visualization ========

figure('position', [1 1 1200 350]);
plot(t_sol, exp(u_sol));
ylim([-0.2 1.2]);
xlabel('time');
ylabel('y(t)');

% End of script


