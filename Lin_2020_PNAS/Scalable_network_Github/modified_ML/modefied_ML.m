
% This code simulates modified May-Leonard model with three nodes. 
% For equation and network diagram, see FigureS2A and Supplementary Method.
% For integration scheme, see 'ode45' function for description.
% It is used to generate the plot in Figure S2B, S3C, S2D.
% WHLin, 2020 Jan

% This simulation is just proof-of-priciple. It is not as accurate 
% as the other one, which apply log-transformed first and hence 
% the simulation closed to 0 is more accurate. 

% ======== Simulation parameters ========

Tstep = 0.1;      % step size of time points 
Tmax = 500;       % maximal time step
optA = odeset('RelTol', 10^(-8),'AbsTol', 10^(-8));  % error tolerance

% Note: the Tstep here is NOT the real integration step size,
% see 'ode45' Matlab document for how integration step size is 
% determined from error tolerance. 

x0 = [2 1 0.5]';    % initial condition
y0 = x0/sum(x0);    % nomalized initial condition 

% ======== Model parameters ========

r = [1.5 1.2 0.9];    % growth parameters
a = 1.1;              % influx coefficients
c = [1 0.8 1.3];      % efflux coefficients

% ======== Model description ========

% Growing trajectory X(t)

N = @(x) x(1)+x(2)+x(3);
y1 = @(x) x(1)/N(x);
y2 = @(x) x(2)/N(x);
y3 = @(x) x(3)/N(x);

F1 = @(x) x(1) * r(1) * ( a - c(1)*y1(x) - c(2)*y2(x) - c(3)*y3(x) );
F2 = @(x) x(2) * r(2) * ( a - c(3)*y1(x) - c(1)*y2(x) - c(2)*y3(x) );
F3 = @(x) x(3) * r(3) * ( a - c(2)*y1(x) - c(3)*y2(x) - c(1)*y3(x) );

mu = @(x) F1(x) + F2(x) + F3(x);
F = @(x) [F1(x) F2(x) F3(x)]';

% Trajectory on rescaled space Y(t)
G = @(y) F(y) - mu(y)*y;

% Integration the solution
Ft = @(t,x) F(x);
[t_sol, x_sol] = ode45(Ft,[0:Tstep:Tmax], x0, optA);

Gt = @(t,y) G(y);
[t_sol2, y_sol] = ode45(Gt,[0:Tstep:Tmax], y0, optA);

% ======== Data visualization ========

figure;
subplot(211);
semilogy(t_sol, x_sol);

subplot(212);
plot(t_sol2, y_sol);

figure;
plot3(y_sol(:,1), y_sol(:,2), y_sol(:,3));
xlabel('Y1');  ylabel('Y2'); zlabel('Y3');


