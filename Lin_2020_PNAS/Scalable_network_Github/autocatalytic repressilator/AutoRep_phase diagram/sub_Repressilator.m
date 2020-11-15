
% This code simulates autocatalytic represillator with four nodes. 
% For integration scheme, see 'ode45' function for description.

% This code is modified for the purpose of simulation across various
% parameter regime for theta and K.

% Input parameters: 
% - yIni: initial condition
% - ParaSet: parameter set for the autocatalytic repressilator 
%   ParaSet(1) = b2;
%   ParaSet(2) = b3;
%   ParaSet(3) = b4;
%   ParaSet(4) = csyn;
%   ParaSet(5) = d2;
%   ParaSet(6) = d3;
%   ParaSet(7) = d4;
% - theta: Hill coefficient parameter 
% - K: repression strength parameter
% - simuT: simulation time
% - StepSize: step size of the solution trajectory

% Note: the step size here is NOT the real integration step size,
% see 'ode45' Matlab document for how integration step size is 
% determined from error tolerance. 

% Ouput result: 
% - LTGR: long-term growth rate, calculated from the second half of y(t)
% - y_avg: y-trajectory for averaging, obtained from the second half of y(t)


function [LTGR, y_avg] = sub_Repressilator(yIni, ParaSet, theta, K, simuT, StepSize)

% ======== Simulation parameters ========

Tstep = StepSize;       % step size of time points
Tmax = simuT/StepSize;  % maximal time step
optA = odeset('RelTol', 10^(-4),'AbsTol', 10^(-4));  % error tolerance

y0 = yIni;            % Initial condition
y0 = y0/sum(y0);      % Normalized initial condition

% ======== Model parameters ========

b = [NaN ParaSet(1:3)];
csyn = ParaSet(4);
d = [NaN ParaSet(5:7)];

% ======== Model description ========

% See Supplementary Materials for flux fucntion description.
% Below flux functions are written for rescaled system, representing J(Y). 
% Note that N=1 in rescaled system and is omitted in the code.

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

time_step_rg = floor(Tmax*[0.5 1]);  % Range of time step for averaging long-term growth rate

t_avg = t_sol(time_step_rg(1):time_step_rg(2));
y_avg = y_sol(time_step_rg(1):time_step_rg(2),:);

% Obtain instantaneous growth rate mu(t)

mu_vec = NaN(size(t_avg,1),1);  

for j = 1:size(t_avg,1)
    mu_vec(j) = mu( y_avg(j,:) );
end

% Calculate long-term growth rate
LTGR = mean(mu_vec);


% Visualizing solution trajectory 
% t_plot = t_avg;
% y_plot = y_avg;

% figure;
% subplot(131); plot(t_plot, y_plot, '.-');
% subplot(132); plot3(y_plot(:,2), y_plot(:,3),  y_plot(:,4), '.'); 
% subplot(133); plot(t_plot, mu_vec, '.-');

% End of the code

