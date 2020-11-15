
% This script analyze simulated dataset and investigate how 
% optimal long-term growth rate can be predicted by energy balance between
% catabolic flux and biosynthesis flux. 

% There are three datasets (in folders PV10, PV20, PV30). These dataset
% were simulated across the same parameter regime; the are only differ
% in the ATP production value (PV) for catabolism. 

% PV is the number of ATP generated when catobolise one amino acid 
% (m7 to m10). For example, PV=20 indicates catabolizing one amino acid
% generates 20 ATP. Intuitively, the higher the PV, the more efficient 
% for the cell to generate energy using catabolism. 

% ======== Calculating long-term growth rate ========

% The simulations are performed by varying two parameters: 
% (1) k_synthesis rate (k_syn)
% (2) k_catabolic rate (k_cata) 
% The parameter k_syn is the rate for polymer synthesis reaction
% The parameter k_cata is the rate for catabolizing amino acid for energy
 
% By varying both parameters, we have a 2-dimensional phase plot 
% for a 2-dimensional parameter regime (see Figure 3E). The anlaysis is
% perform on the results of the phase plot. 

% For a given fixed k_syn, there is an optimal value of k_cata
% which gives highest long-term growth rate. This value can be seen from
% the color map in Figure 3E, (left panel). The highest long-term growth rates 
% for are externel nutrient levels are plotted as a continuous white curve. 
% This is an empirical law between growth rate and k_cata. 

% ======== Studying catabolic fraction ===========

% To see how this empirical law can be associated with energy balance, 
% we focus on the energy balance (characterized by catabolic flux fraction)
% in the same regime. Catabolic flux fraction (defined in main text) is the
% fraction of amino acid that is utilized for generating ATP. 

% Heuristically, the optimal number for catabolic fraction depends on the 
% energy balance, and is related with the ATP production stoichiometry
% (PV). For this specific biosynthesis model, the formula is 
% optimal catabolic fraction = 5 / (5 + PV - 2)
% (see main text for explanation). 

% When we simulate the biosynthesis mode, we obtain the information
% not only for long-term growth rate but also for fluxes. This enable us to
% correlate optimal long-term growth rate with optimal cataolic flux
% fraction. 

% To see if this optimal catabolic fraction really coincidence with optimal
% long-term growth rate, we plot the curve where the cell reaches this number 
% in the same regime (red curve of Figure 3E, right panel). This curve is
% largely coincidence with the white curve in Figure 3E, left panel. 

% The analysis mentioned above are used in Figure 3E and Figure S4.

% =====================================================================

% (1) Set ATP production value (For the dataset here, choose either 10, 20 or 30)
production_value = 30;

cat_gain = production_value - 2;   % net catabolic gain 
syn_cost = 5;                     % net synthesis cost
expected_frac = syn_cost/(cat_gain + syn_cost);   % expected optimal catabolic fraction

% =====================================================================

% (2) Loading the simulated dataset 
data_folder = strcat('PV', num2str(production_value) ,'/');
data_name = 'dataA.mat';
temp = load(strcat(data_folder, data_name));
temp2 = temp.data;

LTGR = temp2.LTGR;   % The 2D phase plot of long-term growth rate 
Js = temp2.Js;       % Equilibrium flux fraction for all simulated conditions

% =====================================================================

% (3) Analyze the phase diagram of long-term growth rate

n1 = size(LTGR, 1);  % the number of different \theta_R in simulation
n2 = size(LTGR, 2);  % the number of different nutrient levels in simulation

opt_rec = NaN(n2,3);  % optimal k_syn record 

% For each k_catabolic rate, find optimal k_syn which gives highest growth
% rate. Each row correspond to one k_cata value
% Column [1]: value of maximal long-term growth rate
% Column [2]: index of k_syn where maximal long-term growth rate is achieved.
% Column [3]: index of k_cata
% Column [4]: smoothed location for the optimal k_syn

for cata_ind = 1:n2
    
    [val,ind] = max(LTGR(:,cata_ind));
    opt_rec(cata_ind,1:3) = [val ind cata_ind];
    
end

% Fit the data with polynomial curve
fit_coef = polyfit(opt_rec(:,3), opt_rec(:,2), 6);
opt_rec(:,4) = polyval(fit_coef, opt_rec(:,3));

% =====================================================================

% (4) Analyze the phase diagram of catabolic fluc fraction

polyN = 120;   % number of monomers for one polymer
catJ = sum(Js(:,:,12:15),3);    % Catabolic flux (relative magnitude)
synJ = polyN*Js(:,:,16);        % Synthesis flux (relative magnitude)
catJf = catJ./(catJ+synJ);      % Catabolic flux fraction (obtaining from simulation)

% We define OBJ = |catabolic flux fraction - expected optimal flux fraction| 
% as an objective fucntion. That is, when OBJ=0 the system reaches 
% optimal flux fraction. 

opt_recB = NaN(n2,3);   % optimal catabolic flux fraction record

% For each k_cata, find optimal k_syn which gives highest growth
% rate. Each row correspond to one nutrient level
% Column [1]: value of that minimize OBJ
% Column [2]: index of k_syn where minimal OBJ is achieved
% Column [3]: index of k_cata.
% Column [4]: smoothed location for the k_syn where OBJ is minimized

for cata_ind = 1:n2
    
    objective_fun = abs( catJf(:,cata_ind) - expected_frac ) ;
    
    [val,ind] = min(objective_fun);
    opt_recB(cata_ind,1:3) = [val ind cata_ind];
  
end

% Fit the data with polynomial curve
fit_coefB = polyfit(opt_recB(:,3), opt_recB(:,2), 6);
opt_recB(:,4) = polyval(fit_coefB, opt_rec(:,3));

% =====================================================================
% (5) Data visualization 

% Generate a customized colormap 'cmap1'

cmap1 = jet;
cmap1(1:15,3) = linspace(0,1,15);
cmap1(1:10,1) = linspace(0,0.3,10);
cmap1(10:25,1) = linspace(0.3,0,16);

% Plot figures 

figure('position', [1 1 400 180]);

xticks = [1 29 59];
xticklabels = [0.1 1 10];
yticks = [1 31 61];
yticklabels = [1 10 100];

%

subplot(121);
imagesc(LTGR); hold on;

set(gca,'YDir','normal')
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);

colormap(cmap1);
contour(LTGR, 'ShowText', 'on', 'color', 'k'); hold on;

plot(opt_rec(:,3), opt_rec(:,4), '-', 'color', 'w', 'LineWidth', 1.5);  hold off;

xlabel('k_{catabolism}');
ylabel('k_{polymer synthesis}');


%
subplot(122);

contour(catJf, 'ShowText', 'on', 'color', 'k'); hold on;
set(gca,'YDir','normal');
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);

plot(opt_recB(:,3), opt_recB(:,4), '-', 'color', 'r');  hold off;

xlabel('k_{catabolism}');
ylabel('k_{polymer synthesis}');

% End of the script

