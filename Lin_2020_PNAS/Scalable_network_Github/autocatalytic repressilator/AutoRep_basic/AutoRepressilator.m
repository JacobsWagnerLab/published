
% This code simulates autocatalytic represillator with four nodes. 
% For equation and network diagram, see Figure2A and Supplementary Method.
% For integration scheme, see 'ode45' function for description.
% WHLin, 2020 Jan

% ======== Simulation parameters ========

Tstep = 0.1;    % step size of time points
Tmax = 1000;    % maximal time step
optA = odeset('RelTol', 10^(-8), 'AbsTol', 10^(-8));  % error tolerance

% Note: the Tstep here is NOT the real integration step size,
% see 'ode45' Matlab document for how integration step size is 
% determined from error tolerance. 

y0 = [0.1 0.3 0.3 0.3];   % Initial condition
y0 = y0/sum(y0);          % Normalized initial condition

% ======== Model parameters ========

b = [NaN 0.5 0.4 0.3];
csyn = 20;
d = [NaN 0.25 0.3 0.2];

theta = 4;  % Hill coefficient
K = 500;    % repression strength

% ======== Model description ========

% See Supplementary Materials for flux fucntion description.
% Below flux functions are written for rescaled system, representing J(Y). 
% Note that N=1 in rescaled system and is omitted in the code.

Jin = @(y) b(2)*y(2) + b(3)*y(3) + b(4)*y(4);
J1 = @(y) Jin(y);

J2 = @(y) csyn* y(1) / (1 + K* ( y(3)^theta ) ) ;
J3 = @(y) csyn* y(1) / (1 + K* ( y(4)^theta ) ) ;
J4 = @(y) csyn* y(1) / (1 + K* ( y(2)^theta ) ) ;

Jd2 = @(y) d(2)*y(2);
Jd3 = @(y) d(3)*y(3);
Jd4 = @(y) d(4)*y(4);

F1 = @(y) J1(y) - J2(y) - J3(y) - J4(y);
F2 = @(y) J2(y) - Jd2(y);
F3 = @(y) J3(y) - Jd3(y);
F4 = @(y) J4(y) - Jd4(y);

F = @(y) [F1(y) F2(y) F3(y) F4(y)]';
mu = @(y) sum(F(y));  % intantaneous growth rate

% The ODE equation on the simplex space
G = @(t,y) F(y) - mu(y).*y;

% Integrating the model
[t_sol, y_sol] = ode45(G,[0:Tstep:Tmax], y0, optA);


% ======== Calculate long-term growth rate ========

% Long-term growth rate is calculated by averaging the 
% instantaneous growth rate, after the rescaled system soltuion Y(t) 
% reach the attractor. Practically, simulate a sufficeint long
% time range and take the average in the end of the simualtion

time_step_rg = floor((Tmax/Tstep)*[0.5 1]);  % Range of time step for averaging long-term growth rate

t_avg = t_sol(time_step_rg(1):time_step_rg(2));
y_avg = y_sol(time_step_rg(1):time_step_rg(2),:);

% Obtain instantaneous growth rate mu(t)

mu_vec = NaN(size(t_avg,1),1);  

for j = 1:size(t_avg,1)
    mu_vec(j) = mu( y_avg(j,:) );
end

% Calculate long-term growth rate
LTGR = mean(mu_vec);


% =========  Visualizing solution trajectory =======
t_plot = t_avg;
y_plot = y_avg;

figure;
subplot(131); plot(t_plot, y_plot, '.-');
subplot(132); plot3(y_plot(:,2), y_plot(:,3),  y_plot(:,4), '.'); 
subplot(133); plot(t_plot, mu_vec, '.-');

% End of the code

