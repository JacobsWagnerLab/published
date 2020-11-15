
% This is the biosynthesis network toy model, contains 24 nodes and 40 fluxes.
% The model is described in Supplementary Materials. The parameter is
% described in Table S3.

% =========== Simulation parameters ===============

Tmax0 = 10;        % maximal time step
Tstep0 = 0.01;     % simu lation time interval
LocAcc = 10^(-6);  % local accuracy
GloAcc = 10^(-6);  % global accuracy
plotflag = 1;      % set plot_flag = 1 for graphic output

sim_condition = [Tmax0, Tstep0, LocAcc, GloAcc, plotflag];

% =========== Model initial condition ===============

% Define 'yini' as the initial condition for simulation

yini_ratioA = [0.4 0.4 0.4 0.14 0.14 0.14 0.5 0.55 0.55 0.55 0.45 0.1];
yini_ratioB = [2 2 2 5 5 5 5 5 6.5 16.5 35 8.5];

yini = [yini_ratioA yini_ratioB];
yini = yini/sum(yini);   

% =========== Model parameters ===============

% External resource concentration
EnvX = 10*[1 1 1];    

% Define proteome partition based on P, Q, R fraction

P_frac = 2;      % Transporter fraction
R_frac = 35;     % Ribosomal fraction
Q_ratio = [3 3 3 3 3 4 10 5];   % House-keeping component ratio
Q_frac = (Q_ratio/sum(Q_ratio))*(100 - (3*P_frac) - R_frac);  % House-keeping component fraction

phi = NaN(1,12);
phi(1:3) = P_frac*[1 1 1];
phi(4:10) = Q_frac(1:7);
phi(11) = R_frac;
phi(12) = Q_frac(8);

phi = phi/sum(phi);

% Note: the remaining parameters are defineed in other subfunctions
% (see define_param.m)


% =========== Simulation start  ===============

% 1. Define flux network topology and flux functions
[J, param_type] = define_fluxes();
        
% 2. Define the stoichiometry matrix S
[S] = define_Smat(phi);

% 3. Define parameters for flux functions
[param_value] = define_param(param_type, EnvX);

% Call subfunction for integration the ODE.

[LTGR, t_sol, y_sol, J_sol] = sub_biosynthesis_model(J, S, param_value, phi, yini, sim_condition);

% Input: 
% - J: the mathematical formula of flux functions
% - S: the stoichiometry matrix
% - param_value: the parameter values
% - phi: proteme partition parameter
% - yini: initial condition
% - sim_condition: simulation condition
% - plotflag: choose 1 for generating plot, choose 0 for not.  

% Ouput: 
% - LTGR: long-term growth rate
% - t_sol: time points of the solution trajectory y(t)
% - y_sol: solution trajectory y(t)
% - J_sol: solution trajecroty of flux J(y(t))
        
LTGR

% End of script

