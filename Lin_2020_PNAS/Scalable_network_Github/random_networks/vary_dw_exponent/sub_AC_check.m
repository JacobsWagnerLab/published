
% This script import a flux set, and remove fluxes step by step until 
% (1) reach a sub-autocatalytic circuit, or 
% (2) the flux set is empty
% See Supplementary Text, Algorithm 6.3 for detail. 
%
% Input: (reaction network for checking the autocatalytic circuit)
% - node_num: number of nodes in reaction network
% - Kini:     flux set of the reaction network
%
% Output:
% - type: type =  1   if an autocatalytic circuit is found. 
%         type = (-1) if no autocatalytic circuit is found.
%
% - AC_size : Number of fluxes in the autocatalytic circuit 
%             If no autocatalytic circuit is founded, AC_size = 0.
%
% - K_auto:  The reduced reaction network which is autocatalytic.
%
%========================================================================

function [type, AC_size, K_auto] = sub_AC_check(node_num, Kini)

flux_num = size(Kini,2);        % number of fluxes in the reaction network
max_reduction_step = flux_num;  % maximal step for iteration of reduction

type = NaN;     % autocatalytic type of the reaction network
K_temp = Kini;  % temporally flux set
K_auto = {};    % the autocatalytic subnetwork (if exist)

% The loop in below performs the search for autocatalytic circuit (AC)
% For each loop, go through all reactions {J_a} and check if 
% mt(phi_a) belongs to dw(K_temp). 

% If for every reaction, mt(phi_a) belongs to dw(K_temp), an AC is found. 
% Otherwise, remove the reactions where dw(J_a) is not in mt(K_temp).
% This resulted in a reduced reaction set 

% Repeat same procedure to the next loop, until an AC is found or 
% the reaction set K_temp is empty. 
% (see Supplementary Text, Algorithm 6.3 for detail) 

for h = 1:max_reduction_step
    
    if (flux_num > 0)
        
        % Call subscript 'reduction.m' to examine the mt(K_temp) and dw(K_temp)
        [Kr, type] = reduction(K_temp, node_num);  
        
        if (type == 1) % found an AC
            
            K_auto = K_temp; 
            break;
            
        elseif (type == (-1)) % reach an empty set
            
            K_auto = {};
            break;
            
        else  % inconclusive, go to the next loop.   
            
            K_temp = Kr;
                                    
        end
        
    else
        
        break;
        
    end    

end

% Report whether the A.C is found or not. 

if (type == 1)
    
    AC_size = size(K_temp,2);

elseif (type == (-1))
    
    AC_size = 0;

else
    
    AC_size = NaN;
    
end

% End of the script s

