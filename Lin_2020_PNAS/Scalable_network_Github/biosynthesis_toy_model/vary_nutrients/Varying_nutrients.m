
% This script loads simulation resulrs (as datasets located in the folder)
% and performs the analysis as in Figure 3B. The dataset is simulated for 
% biosynthesis model with different external nutrient levels (see main text for description) 
% and used to investigate how long-term growth rate is affected by external nutrients. 

% ============================ Load datasets ================================ %

% The dataset is obtained by simulating the biosynthesis model for various nutrients. 

L = 36;   % number of different nutrient levels
n = 3;    % number of nutrient type (s1,s2,s3)

LTGR_vec = zeros(L,n);      % Data vector of long-term growth rate
lim_aa_vec = zeros(L,n);    % Data vector of limiting metabolic building block

% Loading dataset. 

parent_folder = pwd;
folder_name = {'VaryN1', 'VaryN2', 'VaryN3'};

for j = 1:n

    cd(parent_folder);
    cd(folder_name{j});
    temp = load('dataA.mat');
    data = temp.data;
    LTGR_temp = data.LTGR;
    
    % Long-term growth rate of each nutrient conditions 
    LTGR_vec(:,j) = LTGR_temp;
    
    % Equilibrium point of Y(t) in the end of simulation
    ys = data.ys;
        
    % Calculate limiting amino acid concentration:
    % i.e. the lowerst level among four amino acids (m7 to m10)
    lim_aa = min(ys(:,7:10)');
    lim_aa_vec(:,j) = lim_aa';  
        
end


% ==================== Fitting Sigmoidal curve ======================= %

% The code below fit the long-term growth rate data 
% to sigmoidal function (see formula in below)

% Sigmoidal function, with three parameter beta1, beta2, beta3.
% The Hill coeffcient is beta3.
sigfun = @(beta,x) beta(1)* (x.^beta(3)) ./ ((beta(2)^beta(3)) + x.^beta(3));

beta_ini = [2 1 1];          % Initial condition for fitting parameter 'beta'
Lscale = logspace(-2,2,L)';  % range of nutrient levels

LTGR_fit = zeros(L,n);   % Sigmoidal curve fit to long-term growth rate curve 
                         % for different nutrient conditions

for j = 1:n
    
    beta_j = nlinfit(Lscale, LTGR_vec(:,j), sigfun, beta_ini);    
    LTGR_fit(:,j) = sigfun(beta_j, Lscale);
    
end

% ======================== Data visualization ========================= %

% The figures in below are used in Figuse 3B.

figure('position', [1 1 700 300]);

subplot(121);

xaxis = logspace(-2,3,36);
semilogx(xaxis, LTGR_vec(:,1), 'or', xaxis, LTGR_vec(:,2), 'ob', xaxis, LTGR_vec(:,3), 'og');     
legend('s_1','s_2', 's_3');
hold on;    

semilogx(xaxis, LTGR_fit(:,1), '-r', xaxis, LTGR_fit(:,2), '-b', xaxis, LTGR_fit(:,3), '-g'); 
hold off;

xlim([0.01 1000]);
xlabel('Nutrient level');
ylabel('Long-term growth rate');

%

subplot(122);

semilogx(xaxis, lim_aa_vec(:,1), 'or', xaxis, lim_aa_vec(:,2), 'ob', xaxis, lim_aa_vec(:,3), 'og');
legend('s_1','s_2', 's_3');

xlim([0.01 1000]);
xlabel('Nutrient level');
ylabel('Limiting amino acid level (%)');


