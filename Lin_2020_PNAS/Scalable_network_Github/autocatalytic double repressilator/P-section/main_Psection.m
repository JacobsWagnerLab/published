
% This script collects the Poicare section from the autocatalytic 
% doubple repressilator model, for various parameter range of alpha. 
% The result is used in Figure S3A. 

% ======== Simulation parameters for Poincare section ========

alpha_vec = [50:50:500]';
M = size(alpha_vec,1);

% ======== Simulation parameters for ODE model ========

Tstep_pre_run = 0.02;  % time step for pre-simulation
Tmax_pre_run = 5000;   % maximal time span for pre-simulation
Tstep = 0.02;          % time step
Tmax = 500;            % maximal time span

noise_level = 0.01;    % scalable noise level

y0 = [10 0.2 0.15 0.15 0.2 0.15 0.15];  % initial condition
y0 = y0/sum(y0);             % normalized initial condition   

% ======== Simulation start ========

Psec = {};  % saved data points for the poincare section

% Each data struct in 'Psec' is an k-by-2 array. 
% first column: value of alpha (constant)
% second column: projection of points on Poincare section.

parfor j = 1:M
    
    param = alpha_vec(j);
    
    % (1) Pre-run the trajectory to let system converge to attractor
    
    [t_pre, y_pre] = sub_DoubleRep(y0, Tstep_pre_run, Tmax_pre_run, param, noise_level);
    
    % (2) Run long-simulation to obtain attractor
    
    y_pre_run = y_pre(end,:);        
    [t_eq, y_eq] = sub_DoubleRep(y_pre_run, Tstep, Tmax, param, noise_level);
    
    % Use the trajectory Y(t) to find Poincare section
    % (the section number and section position are specified in the subfunction)
    
    Psec{j}.data = sub_Poincare_section(param, y_eq);
    
end

% ======== Data Visualization ========

figure;

for m = 1:M
    
    data_temp = Psec{m}.data;
    
    if ( size(data_temp,1) >1  )
        
        plot(data_temp(:,1), data_temp(:,2), 'k.', 'MarkerSize', 2);
        hold on;

    end
    
end

hold off;

xlabel('\alpha');

% End of the script

