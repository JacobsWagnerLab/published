
% This script defined the value of parameters (total 78 values) in the 
% biosynthesis model. The parameters are ordered by the order of flux function 
% they are used for. The parameters are also listed in Table S3.

% Input: 
% - input: The 78-by-3 matrix specifying the parameter properties. 
%          This is an output from subscript 'define_fluxes.m'.

% - EnvX:  A 1-by-3 vector specifying the environmental nutrient level.
%          For example, EnvX = 10*[1 1 1];

% Ouput: 
% - output: The 78-by-1 array containing the values of 78 parameters.
%           The orders are consistent with the 'input' parameter property.

% ====================================================================== %

function [output] = define_param(input, EnvX)

p = {};  
% 'p' is a cell array that storage parameter or each flux seperately.
% It contains 40 data vectors, each data vector is the parameters set
% of one flux. 

% The following three parameters are used for later definitions. 

KA = 0.005;   % general affinity of ATP
rcat0 = 6;    % general speed for catabolic reaction 
kq0 = 10;     % general kinetic speed for polymerization

%

% In below, parameters are defined for each flux. 
% See Supplementary Mateerials and the subfunction 'define_fluxes.m' 
% for the name and description of these fluxes. Referred to the 
% original flux fucntion TYPE to see where the parameter appears 
% in the flux fucntion.

%============================================================%
 
% (1) Parameter for flux Ain-Cin (TYPE MS1 = 1)
% These are transporter fluxes 

rX = 60*[1 1 1];   % binding rates for transporter
KX = 1*[1 1 1];    % Michaelis-Menten constants for transporter
  
r_import = rX.* (EnvX./ (KX+EnvX));   % Michaelis-Menten constant, k_cat
KA_import = KA;                       % Michaelis-Menten constant, K

p{1} = [r_import(1) KA_import];
p{2} = [r_import(2) KA_import];
p{3} = [r_import(3) KA_import];

%============================================================%

% (2) Parameter for flux S1-S3 (TYPE MS3 = 3)
% These are binary synthesis fluxes for small metabolites.

rS = 120*[1 1 1];    % Michaelis-Menten constants, k_cat
KAS = KA*[1 1 1];    % Michaelis-Menten constants, K
rhoS1 = 0.01;        % ratio parameter between rate constants
rhoS2 = 0.01;        % ratio parameter between rate constants

p{4} = [rS(1) KAS(1) rhoS1 rhoS2];
p{5} = [rS(2) KAS(2) rhoS1 rhoS2];
p{6} = [rS(3) KAS(3) rhoS1 rhoS2];

%============================================================%

% (3a) Parameter for R1F, R2F (TYPE MS3 = 3)
% These are reversible forward-conversion fluxes
 
rF = 250*[1 1];    % Michaelis-Menten constants, k_cat
KAF = KA*[1 1];    % Michaelis-Menten constants, K
rhoF1 = 0.005;     % ratio parameter between rate constants
rhoF2 = 0.005;     % ratio parameter between rate constants

p{7} = [rF(1) rhoF1 rhoF2 KAF(1)];
p{8} = [rF(2) rhoF1 rhoF2 KAF(2)];

%============================================================%

% (3b) Parameter for R1R, R2R (TYPE MS2 = 2)
% These are reversible reverse-conversion fluxes 

rR = 1*[1 1];       % Michaelis-Menten constants, k_cat
KR = 0.04*[1 1];    % Michaelis-Menten constants, K
rhoR = 1*[1 1];     % ratio parameter between rate constants

p{9} = [rR(1) KR(1) rhoR(1)];
p{10} = [rR(2) KR(2) rhoR(2)];

%============================================================%

% (4a) Parameter for SADP (TYPE MS1 = 1) 
% This is the de novo ADP-synthesis flux

rADPsyn = 1;       % Michaelis-Menten constant, k_cat
KADPsyn = 0.05;    % Michaelis-Menten constant, K

p{11} = [rADPsyn KADPsyn];

%============================================================%

% (4b) Parameter for C1-C4 (TYPE MS2 = 2) 
% These are catabolic fluxes that recharges ADP to generate ATP. 
% through consumption of small metabolites.

rcat = rcat0*[1 1 1 1];    % Michaelis-Menten constants, k_cat
Kcat = 0.001*[1 1 1 1];    % Michaelis-Menten constants, K
rhoC = 0.05;

p{12} = [rcat(1) Kcat(1) rhoC];
p{13} = [rcat(2) Kcat(2) rhoC];
p{14} = [rcat(3) Kcat(3) rhoC];
p{15} = [rcat(4) Kcat(4) rhoC];

%============================================================%

% (5) Parameter for polymerization reaction (TYPE Poly = 4) 
% There are 12 polymers and hence 12 synthesis fluxes (Jp1~Jp3, Jq1~Jq8, JR). 
% These polymer synthesis fluxes are proportional through proteome partition ratio,  
% and hence is desbribed by a single synthesis flux.

KZ_poly = 0.002;             % binding equilibrium constant
kq_poly = kq0;               % Universal pre-factor rate
Ka_poly = KA*0.5;            % Michaelis-Menten constants, K
rho_poly = 0.005;            % Universal Sb rho
nb_poly = [30 30 30 30];     % Universal unit of monomer per type

p{16} = [KZ_poly kq_poly Ka_poly rho_poly nb_poly];

%============================================================%

% (6a) Parameter for polymer degradation (type CD = 5) 

DU_poly = 0.005;  % linear degradation rate

for m = 1:12
    
    Jind = m + 16;  % start from p{17}
    p{Jind} = DU_poly;
    
end

%============================================================%

% (6b) Parameter for metabolite degradation (type CD = 5) 

DU_m = 0.005;     % linear degradation rate

for m = 1:12
    
    Jind = m + 28;  % start from p{29}
    p{Jind} = DU_m;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Go through each flux and assign proper parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

row_num = size(input,1); 
output = NaN(row_num, 1);

for r = 1:row_num
    
    tagJ = input(r,1);   % Index of the flux function   
    tagP = input(r,2);   % The order of this parameter in the associated flux function
    temp = p{tagJ};      % The value of the r'th parameter 
    output(r) = temp(tagP);   % write this value in r'th entry of the output array 
    
end

% End of the script



