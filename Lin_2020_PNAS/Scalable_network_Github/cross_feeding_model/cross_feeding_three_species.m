
% This code simulates cross-feeding model with three species. 
% For equation and network diagram, see Figure4 and Supplementary Method.
% For integration scheme, see 'ode45' function for description.
% It is used to generate the plot in Figure S2B, S3C, S2D.
% WHLin, 2020 Jun

Tstep = 1;        
Tmax = 2000;

Km = 0.1;         % Km parameter
r = [1 1 1];      % nature growth rate

a = [0.536 0.486 0.214; 0.266 0.055 0.763; 0.774 0.654 0.182];  % competition matrix
b = [0.2 0.2 0.2];                                              % metabolite secretion rate
c = [0 0.044 0.098; 0.251 0 0.212; 0.436 0.426 0];              % cross-feeding matrix

yini = [0.25 0.25 0.4 0.02 0.03 0.05]'; % initial condition

% ================================================== %

Jr1 = @(y) r(1)* y(1);
Jr2 = @(y) r(2)* y(2);
Jr3 = @(y) r(3)* y(3);

yv = @(y) [y(1) y(2) y(3)]';

J1 = @(y) Jr1(y) * ( 1 - a(1,:) * yv(y) ) ;
J2 = @(y) Jr2(y) * ( 1 - a(2,:) * yv(y) ) ;
J3 = @(y) Jr3(y) * ( 1 - a(3,:) * yv(y) ) ;

% Secretoty fluxes

J4 = @(y) b(1) * y(1);
J5 = @(y) b(2) * y(2);
J6 = @(y) b(3) * y(3);

% Cross-feeding flux

W12 = @(y) c(1,2) * y(1) * ( y(5)/(Km + y(5)) );
W13 = @(y) c(1,3) * y(1) * ( y(6)/(Km + y(6)) );

W21 = @(y) c(2,1) * y(2) * ( y(4)/(Km + y(4)) );
W23 = @(y) c(2,3) * y(2) * ( y(6)/(Km + y(6)) );

W31 = @(y) c(3,1) * y(3) * ( y(4)/(Km + y(4)) );
W32 = @(y) c(3,2) * y(3) * ( y(5)/(Km + y(5)) );


%

F1 = @(y) J1(y) - J4(y) + W12(y) + W13(y);
F2 = @(y) J2(y) - J5(y) + W21(y) + W23(y);
F3 = @(y) J3(y) - J6(y) + W31(y) + W32(y);

F4 = @(y) J4(y) - W21(y) - W31(y);
F5 = @(y) J5(y) - W12(y) - W32(y);
F6 = @(y) J6(y) - W13(y) - W23(y);

F = @(y) [F1(y) F2(y) F3(y) F4(y) F5(y) F6(y)]';
mu = @(y) sum(F(y));

% Project to the simplex space
G = @(t,y) F(y) - mu(y).*y;

% Normalizing initial condition
yini = yini/sum(yini);

% Numerical integartion
optA = odeset('RelTol', 10^(-5),'AbsTol', 10^(-5));
[t_sol, y_sol] = ode45(G, [0:Tstep:Tmax], yini, optA);


% Data visualization

figure;  
subplot(211);
plot(t_sol, y_sol(:,1), 'b-'); hold on;
plot(t_sol, y_sol(:,2), 'g-'); hold on;
plot(t_sol, y_sol(:,3), 'r-'); hold off;

subplot(212);
plot(t_sol, y_sol(:,4), 'b-.'); hold on;
plot(t_sol, y_sol(:,5), 'g-.'); hold on;
plot(t_sol, y_sol(:,6), 'r-.'); hold off;


