
% This script generates a random reaction network and test if it is
% autocatalytic. The flux function is not specified; instead, we generate 
% the maintenance set ramdonly accoring to probability distribution. 

% We specified the random network as following (see Supplementary Material 
% for detail): each reaction has constant numbers of maintenance and downstream
% set. The nodes are randomly chosen according to a probability
% distribution on these node. (That is, some nodes has higher chance to be
% chosen as mt. set; some other nodes has higher chance to be chosen as
% dw. set). The probability is determined by random simulation. 

% Input:
% - mt_num:   number of maintenance set of each reaction flux
% - dw_num:   number of downstream set of each reaction flux
% - node_num: number of nodes in this network
% - flux_num: number of fluxes in this network
% - exponent_param: parameter for exponent of power-law weighting fucntion

% Ouput: 
% - type: if network is autocatalytic, type = 1.
%         if network is not autocatalytic, type = (-1).

% ======================================================================= 

function [type] = sub_main_repeat(mt_num, dw_num, node_num, flux_num, exponent_param)

K = {};  % reaction network

% The data structure of 'K' is a 1-by-m cell (m = number of fluxes) 
% Each struct K{j} contains two lists
% K{j}.mt is the list of maintenance set
% K{j}.dw is the list of downstream set

% =======================================================================

% Generate cdf for choosing the mt set
node_pdf1 = ones(node_num, 1)/node_num;
node_cdf1 = get_cdf(node_pdf1);

% Generate cdf for dw set
[node_pdf2] = power_law_weighting(node_num, exponent_param);
node_cdf2 = get_cdf(node_pdf2);

% (3) Assign mt(K) and dw(K) based on randomly-generated cdfs.

for j = 1:flux_num

    mt_temp = get_sample(node_cdf1, mt_num);
    dw_temp = get_sample(node_cdf2, dw_num);
    
    K{j}.mt = mt_temp;    
    K{j}.dw = dw_temp; 

end

% (4) Check if this random reaction network is autocatalytic.

[type, Ac_size, K_auto] = sub_AC_check(node_num, K);
% - type: type =  1   if an autocatalytic circuit is found. 
%         type = (-1) if no autocatalytic circuit is found.

% End of the script




