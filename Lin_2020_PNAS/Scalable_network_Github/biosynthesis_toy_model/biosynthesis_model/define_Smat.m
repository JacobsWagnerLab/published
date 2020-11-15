
% This script define the stoichiometry matrix of the biosynthesis reaction
% network. The stoichiometry matrix S is an n-by-m matrix, where n is the
% number of nodes, and m is the number of fluxes. Each column in the matrix
% records the stoichiomatry of one reaction flux, with upstream/downstream
% being positive/negative value. (See Supplementary Materials for detail). 

% This script defines S through each reactions (see below). There are 
% 24 nodes and 40 reactions in this model, and hence S is an 24-by-40 matrix. 

% The only required input parameter is 'phi', the proteome partition vector. 
% In this model, the synthesis flux of 12 polymers (node #13-24) depend on 
% the ratio of 'phi', and hence can be combined into a single reaction. 

% ** Note that for this model, we measure all metabolites and polymers in the
% same unit -- biomass (where metabolite monomer = 1 Unit) and assumes
% mass-conservation in reaction. (except import and degradation of monomer)

% That is, the polymerization reaction of N monomers into 1 polymer is written as 
% {monomer (N grams) -> polymer (N grams)}, not {N monomers -> 1 polymer}. 
% Measuring metabolite/polymer in biomass unit facilitate the comparison between 
% node and flux magnitudes. For differential equation, this formulation is
% equivalent to the traditional method, the difference on coefficient can
% be absorbed into the prefector of flux fucntions. 

% =========================================================================

function [S] = define_Smat(phi)

n_num = 24;         % number of nodes
J_num = 40;         % number of reactions 
S = zeros(n_num, J_num);    % Stoichiometry matrix

% ========================================================================= 

% Following parameters are use for later definitions. 

unit_cost = 3;              % cost of ATP per polymerization step
cat_value = 20*[1 1 1 1];   % ATP production stoichiometry, the number of ATP 
                            % generated for cataboliting an building block (#7-10)

unit_num = [30 30 30 30]';  % number of monomer unit of #7-#10
                            % is requird for synthesize a polymer

% ========================================================================= 

% (1) Transporter fluxes, Ain-Cin (consumes ATP): J# = 1-3

for k = 1:3
    
    S(k,k) = 1;
    S(11,k) = -1;
    S(12,k) = 1;
    
end

% (2) Binary synthesis fluxes, S1-S3:  J# = 4-6

S(1,4) = -0.5;
S(2,4) = -0.5;
S(4,4) = 1;
S(11,4) = -S(4,4);
S(12,4) = S(4,4);

S(1,5) = -0.5;
S(3,5) = -0.5;
S(5,5) = 1;
S(11,5) = -S(5,5);
S(12,5) = S(5,5);

S(2,6) = -0.5;
S(3,6) = -0.5;
S(6,6) = 1;
S(11,6) = -S(6,6);
S(12,6) = S(6,6);

% (3a) Reversible forward-conversion fluxes, R1F, R2F : J# = 7-8

S4a = 0.8;
S5a = 0.4;
S5b = 0.4;
S6b = 0.8;
Spd = 0.6;

S(4,7) = -S4a;
S(5,7) = -S5a;
S(7,7) = Spd;
S(8,7) = Spd;

S(5,8) = -S5b;
S(6,8) = -S6b;
S(9,8) = Spd;
S(10,8) = Spd;

% (3b) reversible reverse-conversion fluxes, R1R, R2R : J# = 9-10

S(4,9) = -S(4,7);
S(5,9) = -S(5,7);
S(7,9) = -S(7,7);
S(8,9) = -S(8,7);

S(5,10) = -S(5,8);
S(6,10) = -S(6,8);
S(9,10) = -S(9,8);
S(10,10) = -S(10,8);

% (4a) de novo ADP-synthesis flux, SADP, C2-C4: J# = 11

S(7,11) = -1;
S(12,11) = 1;

% (4b) catabolic fluxes that recharges ADP to generate ATP, 
%  through consumption of small metabolites, C1-C4.  J# = 12-15

S(7,12) = -1;
S(12,12) = (-1) *cat_value(1);    % catalytic value, as define earlier
S(11,12) = cat_value(1);

S(8,13) = -1;
S(12,13) = (-1) *cat_value(2);
S(11,13) = cat_value(2);

S(9,14) = -1;
S(12,14) = (-1) *cat_value(3);
S(11,14) = cat_value(3);

S(10,15) = -1;
S(12,15) = (-1) *cat_value(4);
S(11,15) = cat_value(4);

% ========================================================================= 

% (5) polymerization reaction flux. This flux combines 12 synthesis fluxes  
% J# = 16

poly_ind = 12;  % index flag for polymer (starts at node #13) 

total_unit = sum(unit_num);   % total number of monomers to synthesize one polymer

S(7, 16) = (-1)* unit_num(1);  
S(8, 16) = (-1)* unit_num(2);    
S(9, 16) = (-1)* unit_num(3);
S(10, 16) = (-1)* unit_num(4);

S(11, 16) = (-1)* unit_cost * total_unit;  % total number of ATP consumption
S(12, 16) = unit_cost * total_unit;        % total number of ADP production

for m = 1:12
    
    frac = phi(m);  % fraction of polymer sysnthesis, between 0 and 1
    
    % Total biomass that becomes the polymer
    % The polymer is not measured by number, but by biomass (using monomer as unit)
    S(poly_ind + m , 16) = frac * total_unit;        

end

% ========================================================================= 

% (6a) degradation fluxes from polymer back to monomer  

polyD_Jind = 16;  % index flag for degradation flux (starts from J# = 17)

for m = 1:12
    
    Jind = polyD_Jind + m;
    Xind = poly_ind + m;
    S(Xind, Jind) = (-1) * total_unit;  
    S(7:10, Jind) = unit_num;
    
end

% (6b) degradation fluxes of monomer (monomer degraded into void)

mD_Jind = 28;  % index flag for degradation flux (starts from J# = 29)

for m = 1:12
    
    mind = m;  
    S(mind, mD_Jind+m) = (-1);
        
end

% =========================================================================

% End of the script


