
% This script generates the phase plot in Figure 2D, E 
% by varying parameters of theta (Hill coefficient)
% and K (repression strength). 

% The phase plot is for observing how growth modality change
% in different parameter regime, and how long-term growth rate 
% is also vary in these regime

% The diameter (Diam) of the limit cycle
% is defined by the maximal distance between any two points 
% on the limit cycle. The diameter of a fixed point is zero. 

%%====== Model parameter ======================================%%

% Common initial condition
xIni = [0.1 0.3 0.3 0.3];

% Common parameter set 
% ParaSet(1) = b2;
% ParaSet(2) = b3;
% ParaSet(3) = b4;
% ParaSet(4) = csyn;
% ParaSet(5) = d2;
% ParaSet(6) = d3;
% ParaSet(7) = d4;

ParaSet = [0.5 0.4 0.3 20 0.25 0.3 0.2];

%%======= Simulation parameter ================================%%


simuT = 1000;     % Simulation time
step_size = 0.1;  % Step size of the reported trajectory 
num_theta = 61;   % number of different theta in the phase diagram
num_K = 61;       % number of diffent K in the phase diagram

% Simulation for parameter of theta between 1 and 4.
theta_vec = linspace(1, 4, num_theta);
% Simulation for parameter of K between 10^1 and 10^4.
K_vec = logspace(1, 4, num_K);

% Data matrices for phase diagram
LTGR_mat = NaN*ones(num_theta, num_K);   % long-term growth rate
Diam_mat = NaN*ones(num_theta, num_K);   % diameter of limit cycle

% Start simulation

parfor j1 = 1 : num_theta
    
    for j2 = 1 : num_K

        theta = theta_vec(j1);
        K = K_vec(j2);
        
        % Simulation of long-term growth rate and obtaining trajectory
        [LTGR, y_avg] = sub_Repressilator(xIni, ParaSet, theta, K, simuT, step_size)        
    
        % Calculating the diameter for limit cycle
        [Diam_LC] = FindDiam(y_avg);
        
        LTGR_mat(j1,j2) = LTGR;
        Diam_mat(j1,j2) = Diam_LC;
        
    end
    
end

%%======= Data Visulization ================================%%

figure;

subplot(121);  imagesc(LTGR_mat); colorbar;  axis off; 
subplot(122);  imagesc(Diam_mat);  colorbar;  axis off; 
colormap(jet);

% End of simulation


