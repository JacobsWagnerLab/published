
% This script apply function TYPE to generate specific flux function.
% This is done by string replacement. 

% For example, the original template for Constant Degradation TYPE (CC)
% is a linear function 'q(1)*u(1)', apply this function for parameter p(i)
% and node y(j) resulted in 'p(i)*y(j)'.

% This script is repetitively called by subfunction 'define_flux.m'. 
% Every time this function is called, it generates a output string
% (which specifies specific flux function)

% Input: 
% - input: The input template string, belongs to one of the TYPE. 
% - yind_new: Specifies the nodes that involved in this flux function
% - flux_type: TYPE of the flux function

% Output:
% - ouput: The output string; specific formula of the flux fucntion.
% - pind_current: current index for parameter

% ========================================================================

function [output, pind_current] = ind_assignV(input, yind_new, flux_type)

% Number of parameters for each TYPE: 
%(i) MS1, (ii) MS2, (iii) MS3, (iv) Poly, (v) CD
param_number = [2,3,4,8,1]; 

% Obtaining slots for new parameters. 
pind_new = GetNewInd( param_number(flux_type) );

% Initial index to be assigned 
yind_ini = (1: length(yind_new));
pind_ini = (1: length(pind_new));

% ========================================================================
% Assign index of variables. These are string replacement operation. 

output = input;

for k = 1: length(yind_ini)
    
    old_k = strcat('u(', num2str(yind_ini(k)), ')');
    new_k = strcat('y(', num2str(yind_new(k)), ')');
    
    output = strrep(output, old_k, new_k);  % replacing for node variable
    
end

for k = 1: length(pind_ini)
    
    old_k = strcat('q(', num2str(pind_ini(k)), ')');
    new_k = strcat('p(', num2str(pind_new(k)), ')');
    
    output = strrep(output, old_k, new_k);  % replacing for parameters
    
end

pind_current = pind_new;  % increase the number of current parameter index

% ========================================================================
% Jcount: a global variable specifiying the index of currnet flux function

% param:  a global variable with k-by-3 format. 
% Each row represent one parameter (associated a flux function)
% column [1]: Index of flux function being associated
% column [2]: The rank (first, second, etc) of this parameter in 
%             the associated flux function.
% column [3]: The type of the associated flux function.

global Jcount;
global param;

for k = pind_new(1:end)
    param(k,1) = Jcount;
    param(k,2) = k - pind_new(1) + 1;
    param(k,3) = flux_type;
end
    
Jcount = Jcount + 1;

% End of the script

