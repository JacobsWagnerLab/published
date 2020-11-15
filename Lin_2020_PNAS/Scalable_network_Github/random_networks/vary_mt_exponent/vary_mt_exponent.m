

% This script investigates how the statistical distribution on 
% downstream set affect the autocatalytic probability for random networks.
% (P_auto). The result is used to generate Figure S6A.

% The phase diagram of P_auto is simulated, by varying two variables 
% (i)  number of reaction fluxes in the system
% (ii) The exponent of power-law function, used in weighting the mt. set
% (see main text for description)
%=======================================================================

% Parameter for random network

mt_num = 3;     % number of maintenance set of each reaction flux
dw_num = 2;     % number of downstream set of each reaction flux
node_num = 100;    % number of nodes in this network
repeatN = 15;      % number of repeats for simulation  

% parameter ranges
exponent_param_vec = [0:0.05:3]';    % range of exponents 
K_vec =[0:4:200]';                   % range of number of reaction fluxes 

n1 = size(exponent_param_vec, 1);
n2 = size(K_vec, 1);

% =====================================================================

P_auto = NaN(n1,n2);

parfor j1 = 1:n1
    j1
    
    for j2 = 1:n2
        
        exponent_param = exponent_param_vec(j1);
        K_num = K_vec(j2);

        [AC_percent] = node_num_statistics(mt_num, dw_num, node_num, K_num, exponent_param, repeatN);
        P_auto(j1,j2) = AC_percent;
        
    end
end

figure;  imagesc(P_auto); hold on;


