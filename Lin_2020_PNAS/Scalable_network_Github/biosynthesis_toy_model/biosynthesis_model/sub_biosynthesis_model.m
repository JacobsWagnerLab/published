
% This script assemble the ODE model (from flux functions) for numerical
% integrating the biosynthesis network model. 

% Input: 
% - J: the mathematical formula of flux functions
% - S: the stoichiometry matrix
% - param_value: the parameter values
% - phi: proteme partition parameter
% - y0: initial condition
% - sim_condition: simulation condition
% - plotflag: choose 1 for generating plot, choose 0 for not.  

% Ouput: 
% - LTGR: long-term growth rate
% - t_sol: time points of the solution trajectory y(t)
% - y_sol: solution trajectory y(t)
% - J_sol: solution trajecroty of flux J(y(t))

% ====================================================================

function [LTGR, t_sol, y_sol, J_sol] = sub_biosynthesis_model(J, S, param_value, phi, y0, sim_condition)

row_num = 24;   % number of nodes
col_num = 40;   % number of 

% ====================================================================
%    Model description 
% ====================================================================

Fvec = [];   % The ODE dx/dt = F(x) = SJ(x)

%%% (1) The loop below multiplies the matrix SJ algebraically 
%       and combine the formula in string format

for row = 1:row_num
    
    temp = '';
    
    for col = 1: col_num
        
        if ( S(row, col) ~= 0 )  % some flux connected
            
            Sja = num2str(S(row, col)); 
            Sja = strcat('(', Sja, ')'); 
            
            Ja = strcat('(', J{col}, ')');            
            temp = strcat(temp, '+', Sja, '*', Ja);    

        end 
        
    end
    
    Fvec{row} = temp;
    
end

% ====================================================================
%  (2) The loop below format the string into F(y,p), 
%      a multivariate multidimensional function 
%      with vector argument y and vector argument p.

%      The argument y will stay to be the variable, 
%      while the argument p will be substitute with parameter values. 


Fvec_str = '';  % A large string containing the formula for F(y,p)

Fvec_str = strcat(Fvec_str, '@(y,p)[');
wspace = {' '};
prime_char = strcat(39);  % ascii code for prime symbol (')

for r = 1:row_num    
    Fvec_str = strcat( Fvec_str, Fvec{r} );
    
    if (r < row_num)
        Fvec_str = strcat( Fvec_str, wspace);
    end
    
end

Fvec_str = strcat(Fvec_str, ']', prime_char);

% Use the built-in fucntion 'str2func' to convert the string
% into the a annomyous Matlab function 

Fyp = str2func(Fvec_str{1});

% =====================================================================

%   (3) In addition, assemble the formula of flux functions into 
%    a large string. This is for substituting the solution Y(t) 
%    back to the flux fucntions and analyzing the flux magnitude. 

Jvec_str = '';   % the large string containing all flux function formula

Jvec_str = strcat(Jvec_str, '@(y,p)[');

for c = 1:col_num    
    Jvec_str = strcat( Jvec_str, J{c} );
    
    if (r < col_num)
        Jvec_str = strcat( Jvec_str, wspace );
    end
    
end

Jvec_str = strcat(Jvec_str, ']', prime_char);

% A multivariate, multidimensional function for all fluxes
Jyp = str2func(Jvec_str{1});   


% =====================================================================
%   Change of variables of ODE 
% =====================================================================

% To improve the accuracy and numeric stability, transform the ODE
% into log space

F = @(y) Fyp(y, param_value);   % original ODE
mu = @(y) sum(F(y));
G = @(y) F(y) - mu(y)*y;        % Project to the simplex space

% Perform log-transform into function of u

w = @(u) exp(u);
GLog = @(u) (1./w(u)) .* G(w(u));
GLogt = @(t,u) GLog(u);


% =====================================================================
%    Simulation condition
% =====================================================================

Tmax = sim_condition(1);     % maximal imulation time span
Tstep = sim_condition(2);    % simulation time interval

% Setting error tolaerance in simulation
optA = odeset('RelTol',sim_condition(3),'AbsTol', sim_condition(4));

% Setting initial condition and transformed in log space
y0 = y0/sum(y0);
u0 = log(y0); 

% Integrate the ODE
[t_sol, u_sol] = ode45(GLogt,[0:Tstep:Tmax], u0, optA);

y_sol = exp(u_sol);   % tranform the ODE from u(t) back to y(t)


% =====================================================================
%    Additional analysis
% =====================================================================

% Further analysis based on the simulation result

% (1) Calculate long-term growth rate

y_size = size(t_sol,1);    % size of solution vector y(t)  
fracI = [0.8 0.99];        % fraction for averaging the long-term growth rate

Rg = floor(y_size*fracI);   % time range for averaging the long-term growth rate

mu_vec = zeros(Rg(2)-Rg(1)+1, 1);

for tj = Rg(1):Rg(2)
    
    mu_vec(tj-Rg(1)+1, 1) = mu(y_sol(tj,:));   %calculate instantaneous growth rate 

end

LTGR = mean(mu_vec);


% (2) Calculate trajectory of flux magnitudes

J_sol = zeros(y_size, col_num);

for ts = 1:y_size        
    J_sol(ts,:) = Jyp(y_sol(ts,:), param_value);
end


% =====================================================================
%   Plot solution trajecotories
% =====================================================================

title_text = {};
title_text2 = {};
title_text3 = {};
title_text4 = {};

title_text{1} = 'metabolites \newline m_1 to m_3';
title_text{2} = 'metabolites \newline m_4 to m_6';
title_text{3} = 'building blocks \newline m_7 to m_{10}';
title_text{4} = 'ATP (blue) \newline ADP (red)';
title_text{5} = 'transporter, enzyme \newline P, Q_1';
title_text{6} = 'catabolic enzymes \newline Q_6, Q_7';
title_text{7} = 'ribosome, RNA \newline R, Q_8';

title_text2{1} = 'transporter influxes';
title_text2{2} = 'primary synthesis';
title_text2{3} = 'reversible fluxes';
title_text2{4} = 'ADP syn. fluxes';
title_text2{5} = 'catabolic fluxes';

title_text3{1} = 'transporter synthesis';
title_text3{2} = 'metabolic \newline EZ synthesis';
title_text3{3} = 'catabolic \newline EZ synthesis';
title_text3{4} = 'ribosome, genetic. \newline synthesis';

title_text4{1} = 'polymer degradation';
title_text4{2} = 'monomer degredation';


%%%%%%=========================================================%%%%%%

if (sim_condition(5) == 1)

% I, I2, I3, I4 are text arrays for the title of the plot
I = {};
I2 = {};
I3 = {};
I4 = {};
    
I{1} = [1:3];
I{2} = [4:6];
I{3} = [7:10];
I{4} = [11:12];
I{5} = [13:20];
I{6} = [21:22];
I{7} = [23:24];

I2{1} = [1:3];
I2{2} = [4:6];
I2{3} = [7:10];
I2{4} = [11];
I2{5} = [12:15];

I3{1} = [1:3];
I3{2} = [4:8];
I3{3} = [9:10];
I3{4} = [11:12];

I4{1} = [17:28];
I4{2} = [29:40];


J_sol_poly = zeros(y_size, 12);

for m = 1:12

    J_sol_poly(:,m) = phi(1,m) * J_sol(:,15);    

end

figure('position', [1, 1, 1600, 800]);

for g = 1:7
    
    subplot(3,7,g); 
    semilogy(t_sol(:,1), y_sol(:,I{g}), 'o-' );  
    title(title_text{g});
end

for g2 = 1:5
        
    subplot(3,7, g2+7); 
    semilogy(t_sol(:,1), J_sol(:,I2{g2}), 'o-' ); 
    title(title_text2{g2});
    
end

for g3 = 1:4
        
    subplot(3,7, g3+12); 
    semilogy(t_sol(:,1), J_sol_poly(:,I3{g3}), 'o-' ); 
    title(title_text3{g3});
    
end

for g4 = 1:2
        
    subplot(3,7, g4+16); 
    semilogy(t_sol(:,1), J_sol(:,I4{g4}), 'o-' ); 
    title(title_text4{g4});
    
end

end

