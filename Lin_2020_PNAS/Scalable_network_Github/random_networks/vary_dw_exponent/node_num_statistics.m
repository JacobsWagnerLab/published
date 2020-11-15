
% This script repetitively generates multiple random networks
% and reports the percentage of these networks are autocatalytic.

% Input:
% - mt_num:   number of maintenance set of each reaction flux
% - dw_num:   number of downstream set of each reaction flux
% - node_num: number of nodes in this network
% - flux_num: number of fluxes in this network
% - exponent_param: parameter for exponent of power-law weighting fucntion
% - repeatN:  number of repeats for simulation  

function [AC_percent] = node_num_statistics(mt_num, dw_num, node_num, K_num, exponent_param, repeatN)

record = NaN(repeatN, 2);

for j = 1:repeatN

    temp = sub_main_repeat(mt_num, dw_num, node_num, K_num, exponent_param);
    record(j,:) = temp;
    
end

% check the percentage of network to be autocatalytic
AC_percent = sum(record(:,1) == 1) / repeatN;


