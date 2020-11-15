
% This script generates weighting factors for N nodes.
% The weighting factor summed up to 1, hence each node could get a 
% different number between 0 and 1. 

% The weighting factor is draw from (truncated) power-law distribution 
% Q(x) = C* x^(-gamma),   x > 0.0001
% with gamma > 0 is the power-law exponent, 
% and C is the normalization constant.

% For gamma = 0, the weighting factor is uniform for all nodes, 
% and hence each node get the same value 1/N.

% For gamma > 0, the weighting factor is not uniform. Some rare nodes 
% will have large weighting factor, while most of the node will have 
% small weighting factors. The distirbution over N nodes follows the power
% law. 

% Practically, in this study the exponent (gamma) is between 0 and 3. 
% We impose a cutoff (for numerical purpose) such that when calculating
% weigting factor no extremely large number will be involved. 
% (see Supplementary Materials for discussion)

% =====================================================================
%
%  Input: 
%  - node_num: number of random variables (iid) to be generated 
%              from power-law distribution
%  - gamma: The exponent parameter for power-law. 
%           In this study, the exponent parmeter is set between 0 and 3. 
%
%  Output:
%  - prob: the weighting facter
%
% =====================================================================
% Example: 
%
% Input:  power_law_pdf2B(10,2)
%
% Output: weighting_factor
%         = [0.0028 0.9361 0.0017 0.0014 0.0026...
%            0.0021 0.0022 0.0078 0.0028 0.0407]';
%
% =====================================================================

function [weighting_factor] = power_law_weighting(node_num, gamma)

% Cutoff for power law (since the power-law formula diverge at 0)
epsilon_cut = 10^(-4);

% Obtaining N random weighting factor draw from power-law distribution
rv = max( rand(node_num, 1), epsilon_cut );

% Power-law distribution 
Q = (-1)* gamma * log(rv);

% un-normalized weighting factor
weighting_factor_temp = exp(Q);

% normalized weighting factor
weighting_factor = weighting_factor_temp /sum(weighting_factor_temp);

