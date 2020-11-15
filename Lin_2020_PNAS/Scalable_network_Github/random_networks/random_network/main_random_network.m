
% This script generates a random reaction network and test if it is
% autocatalytic. The flux function is not specified; instead, we generate 
% the maintenance set ramdonly accoring to probability distribution. 

% We specified the random network as following (see Supplementary Material 
% for detail): each reaction has constant numbers of maintenance and downstream
% set. The nodes are randomly chosen according to a probability
% distribution on these node. (That is, some nodes has higher chance to be
% chosen as mt. set; some other nodes has higher chance to be chosen as
% dw. set). The probability is determined by random simulation. 

% ======================================================================= 

% (1) Define parameter set

node_num = 50;    % number of nodes in this network
flux_num = 70;    % number of fluxes in this network

mt_num = 3;   % Number of maintenance set of each reaction flux
dw_num = 2;   % Number of downstream set of each reaction flux

K = {};   % reaction network

% The data structure of 'K' is a 1-by-m cell (m = number of fluxes) 
% Each struct K{j} contains two lists
% K{j}.mt is the list of maintenance set
% K{j}.dw is the list of downstream set

% ======================================================================= 
% (2) Generate two weighting cumulative distributions functions (cdf)
%     for nodes, which will be used for assigning the maintenance and
%     downstream sets

% Generate cdf for choosing the mt set
[node_pdf1] = power_law_weighting(node_num, 0.5);
node_cdf1 = get_cdf(node_pdf1);

% Generate cdf for choosing the dw set
node_pdf2 = ones(node_num, 1)/node_num;
node_cdf2 = get_cdf(node_pdf2);
    
% ======================================================================= 

% (3) Assign mt(K) and dw(K) based on randomly-generated cdfs.

for j = 1:flux_num

    K{j}.mt = get_sample(node_cdf1, mt_num);
    K{j}.dw = get_sample(node_cdf2, dw_num);
    
end

% (4) Check if this random reaction network is autocatalytic.

[type, AC_size, K_auto] = sub_AC_check(node_num, K);

AC_size






