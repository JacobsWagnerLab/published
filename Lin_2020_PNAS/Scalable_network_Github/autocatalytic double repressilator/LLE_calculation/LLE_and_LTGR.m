
% This is a script for calculating the largest Lyapunov exponent (LLE), 
% see Supplementary Material for detail algorithem. 
% Here, the key parameter being varied is alpha, see Supplementary
% Materials for description of the Autocatalytic double repressilator. 

% Note that calculating LLE involves random perturbation of the trajectory
% and hence requires multiple runs and obtaining statistics.  

% ======== Simulation parameters of the ODE model ========

% The ODE model is described in sub_DoubleRep.m 
% More parameters can be modified in that document.

y0 = [0.21 0.2 0.1 0.1 0.2 0.1 0.1];   % initial condition
y0 = y0/sum(y0);

Tstep0 = 0.1;       % time step for pre-simulation 
Tmax = 200;         % maximal time for pre-simulation
noise_level = 0;    % scalable noise level

% ======== Simulation parameters for calculating LLE ========

param_vec = (100:10:200)';  % Specify the values of parameter (alpha) to be calculated
N = size(param_vec,1);      % Numbers of alpha 
repeatN = 5;                % number of repeated simulation

perturb_scale = 0.1;   % Magnitude of the perturbation
Tstep = 0.01;          % Time step for simulation of Y(t) and purturbed Y(t)
Tstep_LLE = 20;        % Time span for calculating the exponential divergent of LE 
total_step_num = 50;   % Number of perturbation test for each trajectory

% ======== Simulation start ========

% Data matricis with results from different parameters

spec_rec = NaN(N,repeatN);    % the spectrum of LLE for different alpha
GR_rec_ori = NaN(N,repeatN);  % the long-term growth rate for different alpha

parfor w = 1:N
    
    % The entire simulation is repeated for multiple time, 
    % since the simulation involves randomness in choosing perturbation.
    
    for rp = 1:repeatN
        
    param = param_vec(w);   % the parameter alpha
    
    % ===================================================================================   
    % (1) This is a pre-simulation procedure, in order for the trajectory
    %     to be converged to attractor
    
    [t_sol, y_sol, GR_avg] = sub_DoubleRep(y0, Tstep0, Tmax, param, noise_level);   
    
    y_ini = y_sol(end,:);  % Using the ending position Y(t_end) as new initial condition
    
    % ====================================================================================    
    % (2) Calculate LLE for multiple time segment. 
    % For each time segment, we compare unperturbed and perturbed trajectory, 
    % and infer the exponential divergent (Lyapunov exponent)

    % (2-1) Generate unperturbed initial condition 
    y_ori0 = y_ini;
    
    % (2-2) Generate a random perturbation
    perturb_seed = randn(1,7);
    perturb_ini = perturb_scale* perturb_seed / norm(perturb_seed);
       
    % The perturb vector (y_ptb) needs to be normalized with sum(yj)=1
    y_ptb_temp = y_ini + perturb_ini;
    y_ptb0 = y_ptb_temp / sum(y_ptb_temp);
    
    % Temporary data for LE and trajectory    
    Record_LLE = zeros(total_step_num, 3);
    Record_ori = zeros(total_step_num, 7); 
    Record_ptb = zeros(total_step_num, 7);
    
    
    % ==================================================================================== 
    % (3) Simulation for LLE for multiple test
    
    for s = 1 : total_step_num
        
        % (3-1) Simulate two trajectories, one for unpertubed initial condition,
        % the other for a perturbed initial condition
        
        Record_ori(s,:) = y_ori0;
        Record_ptb(s,:) = y_ptb0;
    
        [t_ori, y_ori, GR_ori] = sub_DoubleRep( y_ori0, Tstep, Tstep_LLE, param, noise_level );  
        [t_ptb, y_ptb, GR_ptb] = sub_DoubleRep( y_ptb0, Tstep, Tstep_LLE, param, noise_level ); 
        
        % Calculating the exponential divergence of two trajectory
        
        delta0 = norm(y_ori(1,:) - y_ptb(1,:));
        delta1 = norm(y_ori(end,:) - y_ptb(end,:));
    
        Record_LLE(s,1) = delta0;
        Record_LLE(s,2) = delta1;
        Record_LLE(s,3) = log(delta1/delta0);
        
        % (3-2) Normalize the unpertubed trajectory with sum(yj)=1
        % The normalized vector is used for the initial condition 
        % for the next test.
        y_ori0 = y_ori(end,:);    
        
        % (3-3) Normalized the perturbation, keep the perturbation direction 
        % but reduce the perturbation magnitude for the next test. 
        
        % The simplex space T^(n-1) is a convex space. Therefore, the 
        % new perturbed vector (scaled the delta_y ) remains in T^(n-1)        
        y_ptb_temp = y_ori0 - perturb_scale * (y_ori(end,:) - y_ptb(end,:))/delta1;
        
        % y_ptb_temp should be in T^(n-1). In case of numerical error,
        % normalized by |y_ptb| again to make sure it is in simplex space.
        y_ptb0 = y_ptb_temp / sum(y_ptb_temp);
        
        
    end
    
    % By averaging across multiple perturbation tests, 
    % obtain averaged LLE and averaged long-term growth rate
    
    LLE = mean(Record_LLE(:,3))/Tstep_LLE;   
    LTGR_ori = mean(GR_ori);
    
    % Save into the data matricies
    spec_rec(w,rp) = LLE;
    GR_rec_ori(w,rp) = LTGR_ori;            
    
    end
    
end

% ======== Data Summarization ========

% Spectrum of LLE, across various alpha.
% From multiple simulation, obtaining mean and standard error.
spec_mean = mean(spec_rec');
spec_ste = std(spec_rec')/sqrt(repeatN-1);

% Calculated the long-term growth rate from multilpe simulation.
% This result should be identicle if there is no noise in the system,
% with noise, the LTGR could be slightly different
GR_ori_mean =  mean(GR_rec_ori,2);


% Plotting long-term growth rate and LLE for various alpha

figure;
plot(param_vec, GR_ori_mean,'.-');
xlabel('\alpha');
ylabel('LTGR');

figure;
errorbar(param_vec, spec_mean, spec_ste, 'o');
xlabel('\alpha');
ylabel('largest Lyapunov exponent (LLE)');

% End of the script



