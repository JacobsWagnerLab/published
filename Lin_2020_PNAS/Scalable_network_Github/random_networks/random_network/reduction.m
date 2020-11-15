
% This script performs the search for autocatalytic circuit (AC)
% For each loop, go through all reactions {J_a} and check if 
% mt(phi_a) belongs to dw(K_temp). 

% If for every reaction, mt(phi_a) belongs to dw(K_temp), an AC is found. 
% Otherwise, remove the reactions where dw(J_a) is not in mt(K_temp).
% This resulted in a reduced reaction set 

% Repeat same procedure to the next loop, until an AC is found or 
% the reaction set K_temp is empty. 
% (see Supplementary Text, Algorithm 6.3 for detail) 

% Input: 
% - Kp: reaction network
% - node_num: number of nodes in the reaction network

% Output:
% - Kr: reduce reaction network (which contains an AC)
% - reduction_type: 
%   type =  1   if an autocatalytic circuit is found. 
%   type = (-1) if no autocatalytic circuit is found.

%=======================================================================

function [Kr, reduction_type] = reduction(Kp, node_num)

% (1) scan mt(J) and dw(J) for each flux, 

flux_num = size(Kp,2);

% Generate counting list (0,1,2,... for each node) from maintenance and downstream set
% % For each node, if it belongs to maintenance set of one reaction, add 1. 
% For each node, if it belongs to downstream set of one reaction, add 1. 

% mt_count = zeros(node_num, 1);   % maintenance checking list  
dw_count = zeros(node_num, 1);   % downstream checking list
% 
% % (1a) Go through all reaction, and increase the "node count" of its 
% % maintenance set by 1. 
% 
% for m = 1:flux_num
%     
%     temp = Kp{m}.mt;  % maintenance set for m'th reaction
%     
%     for p = 1:size(temp,1)
%             
%         ind = temp(p);
%         mt_count(ind) = mt_count(ind) + 1;
% 
%     end
%     
% end

% Go through all reaction, and increase the "node count" of its 
% downstream set by 1. 

for m = 1:flux_num
    
    temp = Kp{m}.dw;  % downstream set for m'th reaction
    
    for p = 1:size(temp,1)
    
        ind = temp(p);
        dw_count(ind) = dw_count(ind) + 1;

    end
    
end

% Generate binary checking list of dw(Kp) for nodes. 
dw_rec = (dw_count > 0);

%=======================================================================

% (2) Classify the flux in Kp

% Generate a checking list for each reaction J_a
% If mt(J_a) belongs to dw(Kp), mark 1.   Otherwise, mark 0.
AC_check_list = zeros(flux_num, 1);  


% For each reaction flux, test if all its maintenance set
% are in the downstream set of Kp

for m = 1:flux_num
    
    flux_m = Kp{m}.mt;    % maintenance set of m'th flux 
    flux_m_flag = 1;      % chaecking flag
    
    for p = 1:size(flux_m,1)
    
        temp_node = flux_m(p);
        
        if ( dw_rec(temp_node) == 1 )

            % mt. node found in dw. of Kp
            
        elseif  ( dw_rec(temp_node) == 0 )
            
            % mt. node missing in dw. of Kp
            flux_m_flag = 0;
            
        end
        
    end    
      
    AC_check_list(m) = flux_m_flag;    
    
end

%=======================================================================

% (3) Reprt the reduction type 
%   type =  1   if an autocatalytic circuit is found. 
%   type = (-1) if no autocatalytic circuit is found.

reduction_type = 0;

if ( sum(AC_check_list) == flux_num )

    % Find a sub-AC.    
    reduction_type = 1;
    
elseif ( sum(AC_check_list) == 0 )
    
    % No flux any more.    
    reduction_type = -1;
    
end

%=======================================================================

% (4) Perform the reduction and return 

Kr = {};          % The reduced reaction network
write_flag = 0; 

% Go through each reaction. Remove this reaction if its maintenance set 
% does not belongs to the downstream set of Kp

for j = 1:flux_num
    
    if ( AC_check_list(j) == 1 )
        
        % log into the next flux set
        write_flag = write_flag + 1;
        Kr{write_flag} = Kp{j};        
    
    end
    
end

% End of the script

